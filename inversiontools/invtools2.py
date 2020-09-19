"""Method for identifying temperature inversion layers in a radiosonde sounding.
The core of the method is to fit a piecewise linear function to the profile, then 
estimate the 1-alpha confidence intervals around the slope for each piece. Pieces 
that have a positive slope at the alpha significance level are merged. The method 
is designed to work with GRUAN netcdf radiosonde data, which includes uncertainty
information for each measurement.

GRUAN data processing includes a low-pass filter that reduces the temporal
resolution of temperature from 1s to 10s, but leaves the data in 1s resolution 
format. As the ascent speed is approximately 5 m/s, this means that the effective 
vertical resolution for temperature is 50m. For relative humidity, wind, and altitude,
the vertical resolution is different. (Dirksen et al., 2014). The effective temporal
resolution of the relative humidity depends on temperature (and other things) and as
such is included as its own column. Below the tropopause, the resolution is 10s. Altitude
resolution is 15 seconds, and wind is much lower resolution, at approximately 40s.
Prior to fitting the linear piecewise function, we interpolate to a regular 50 m grid 
starting at ground level. The grid is represented by the parameter dz, and so that 
sensitivity checks can be run is left as a user-specifiable parameter.

Uncertainty in dewpoint temperature: I'm still working on recreating the results of
the paper, but in the meantime, I use the min/max method to get an upper estimate on
the uncertainty. For Td(Ta, RH), the estimate is 
$$U_{T_d} = 1/2*(T_d(T_a+U_{T_a}, RH+U_{RH}) - T_d(T_a-U_{T_a}, RH-U_{RH})$$

compute_inversions:
    applies inversion algorithms to specified netcdf files, including classification.
"""

import xarray as xr
import numpy as np
import pandas as pd
import os
import warnings

warnings.filterwarnings('ignore', 'invalid value encountered in less', RuntimeWarning)
warnings.filterwarnings('ignore', 'invalid value encountered in greater', RuntimeWarning)
warnings.filterwarnings('ignore', 'invalid value encountered in true_divide', RuntimeWarning)
warnings.filterwarnings('ignore', 'divide by zero encountered in true_divide', RuntimeWarning)

def default_params():
    """Default set of parameters. Parameters let the other functions
    know what column names are as well as set options for the classification
    of temperature inversions.
    TODO: Add units, at least in description of default_params
    """
    
    return {'temperature': 'temperature',
           'height': 'altitude',
           'pressure': 'pressure',
           'potential_temperature': None,
           'time': 'date',
           'u_temperature': 'uncertainty_temperature',
           'u_height': 'uncertainty_altitude',
           'dewpoint_temperature': 'dewpoint_temperature',
           'u_dewpoint_temperature': 'uncertainty_dewpoint_temperature',
           'iso_base': False,
           'iso_above': False,
           'iso_within': True,
           'iso_alone': False,
           'na_tol': 0.1,
           'max_embed_depth': 0,
           'min_inv_depth': 75,
           'max_lapse_rate': 0.4,
           'classify_inversions': True,
           'k': 2,
           'dz': 75 # Could set it up so that dz is also potentially a vector
           }

def add_lapse_rate(dataset, params, dewpoint=False):
    """Computes the lapse rate and returns a dataset with lapse rate added
    along with a column with the k=1 uncertainty for the lapse rate. Assumes
    that the altitude or height is provided as the vertical dimension of the
    dataset. Note that here the uncertainty is the uncertainty of the slope of
    the piecewise linear function defined by the sounding data, not the uncertainty
    in the forward difference estimate of the derivative."""
    
    # Note: may need to switch to backward difference to avoid 
    # inconsistency with inversion definitions?
    
    # step one: apply forward difference
    dt = dataset.diff(dim=params['height'], n=1, label='lower').variables[params['temperature']]
    dz = params['dz'] # Here, can alter to use variable spacing if needed. Might be better if
                      # this code is to be used for inversion strength also
    
    u_dt = ((dataset.variables[params['u_temperature']].shift(shifts={params['height']: 1})**2 + 
            dataset.variables[params['u_temperature']]**2)**0.5)[:,0:-1].values

    u_dz = ((dataset.variables[params['u_height']].shift(shifts={params['height']: 1})**2 + 
            dataset.variables[params['u_height']]**2)**0.5)[:,0:-1].values

    n = dataset.coords[params['height']].values.size
    m = dataset.coords[params['time']].values.size

    dt_dz = np.full((m, n), np.nan)
    u_dt_dz = np.full((m, n), np.nan)
    dt_dz[:, :-1] = dt.values/dz
    u_dt_dz[:, :-1] = ((dt/dz)**2 * ((u_dt/dt)**2 + (u_dz/dz)**2))**0.5
    dt_dz = xr.DataArray(dt_dz, dims={params['time']: dataset.coords[params['time']],
                                        params['height']: dataset.coords[params['height']]})
    u_dt_dz = xr.DataArray(u_dt_dz, dims={params['time']: dataset.coords[params['time']],
                                        params['height']: dataset.coords[params['height']]}) 
    if dewpoint:
        dtd = dataset.diff(dim=params['height'], n=1, label='lower').variables[params['dewpoint_temperature']]
        u_dtd = ((dataset.variables[params['u_dewpoint_temperature']].shift(shifts={params['height']: 1})**2 + 
                dataset.variables[params['u_dewpoint_temperature']]**2)**0.5)[:,0:-1].values
        dtd_dz = np.full((m, n), np.nan)
        u_dtd_dz = np.full((m, n), np.nan)
        dtd_dz[:, :-1] = dtd.values/dz
        u_dtd_dz[:, :-1] = ((dt/dz)**2 * ((u_dt/dt)**2 + (u_dz/dz)**2))**0.5
        dtd_dz = xr.DataArray(dtd_dz, dims={params['time']: dataset.coords[params['time']],
                                            params['height']: dataset.coords[params['height']]})
        u_dtd_dz = xr.DataArray(u_dtd_dz, dims={params['time']: dataset.coords[params['time']],
                                            params['height']: dataset.coords[params['height']]}) 
        return dataset.assign({'lapse_rate': dt_dz,
                               'dewpoint_lapse_rate': dtd_dz,
                               'uncertainty_lapse_rate': u_dt_dz,
                               'uncertainty_dewpoint_lapse_rate': u_dtd_dz})
    else:
        return dataset.assign({'lapse_rate': dt_dz,
                           'uncertainty_lapse_rate': u_dt_dz})

    


def find_inversions(data, params):
    """Given a dataframe with a single sounding, find all inversion layers."""

    data.reset_index(inplace=True)
    lr = data['lapse_rate'].values
    ulr = data['uncertainty_lapse_rate'].values
    sign_vector = np.full((len(data),), np.nan)
    # Where lr is 2 sigma above zero, classify as inversion
    # Where lr is 2 sigma below zero, classify as noninversion
    # Where lr - 2 sigma is less than zero but lr + 2 sigma is positive, classify as
    sign_vector[(lr - k*ulr) > 0] = 1
    sign_vector[(lr + k*ulr) < 0] = -1
    sign_vector[((lr - k*ulr) < 0) & ((lr + k*ulr) > 0)] = 0

    def build_layer_df(sign_vector, data):
        """Selects the data at the top and bottom of layers of constant sign
        based on the sign vector, differences the data across the layers, and
        adds columns 'layer_strength' and 'layer_depth'."""
        idx_list = np.atleast_1d(np.argwhere(sign_vector[1:] != sign_vector[:-1]).squeeze())

        idxb = np.array(idx_list[:-1])+1
        idxt = np.array(idx_list[1:])+1
#       # I'll have to figure this part out in a bit,  
#         if len(idx_list) % 2 == 1:
#             idx_list = np.concatenate((idx_list, [len(sign_vector)-1]))

        layer_dict = {'idx_base': idxb,
                      'idx_top': idxt}

        for cc in [cc for cc in data.columns if (cc != 'date')]:
            layer_dict[cc + '_base'] = data.loc[idxb, cc].values
            layer_dict[cc + '_top'] = data.loc[idxt, cc].values

        layer_dict['layer_strength'] = layer_dict[params['temperature'] + '_top'] - layer_dict[params['temperature'] + '_base']
        layer_dict['layer_depth'] = layer_dict[params['height'] + '_top'] - layer_dict[params['height'] + '_base']

        layer_df = pd.DataFrame(layer_dict)
        layer_df = layer_df.dropna(subset=['layer_strength']).reset_index(drop=True)
        if len(layer_df) == 0:
            layer_df = pd.DataFrame(columns = layer_df.columns, index=[0])
        layer_df['date'] = np.array(data.date.values[0])    
        return layer_df

    def merge_iso_layers(sign_vector, params):
        """Based on settings in params, add isothermal
        layers to inversion layers."""

        # Apply rules for including isothermal layers
        if params['iso_alone']:
            sign_vector[sign_vector == 0] = 1
            return sign_vector

        if params['iso_base']:
            iso_idx = (sign_vector[:-1] == 0) & (sign_vector[1:] == 1)
            while np.any(iso_idx):
                sign_vector[:-1][iso_idx] = 1
                iso_idx = (sign_vector[:-1] == 0) & (sign_vector[1:] == 1)

        if params['iso_above']:
            iso_idx = (sign_vector[:-1] == 1) & (sign_vector[1:] == 0)
            while np.any(iso_idx):
                sign_vector[1:][iso_idx] = 1
                iso_idx = (sign_vector[:-1] == 1) & (sign_vector[1:] == 0)

        if not params['iso_alone']:
            sign_vector[sign_vector==0] = -1

        return sign_vector

    def drop_thin_inv_layers(sign_vector, layer_df, params):
        """Inversion layers thinner than the minimum inversion
        depth are dropped. If minimum inversion depth is smaller
        than dz, then this code doesn't do anything."""

        thin_layers = ((layer_df['layer_strength'] > 0) &
                       (layer_df['layer_depth'] < params['min_inv_depth'])).values
        if np.any(thin_layers):
            for x in np.argwhere(thin_layers):
                idxb = int(layer_df.loc[x, 'idx_base']) + 1
                idxt = int(layer_df.loc[x, 'idx_top'])
                sign_vector[idxb:idxt] = -1
        return sign_vector

    def merge_neg_layers(sign_vector, layer_df, params):
        """Thin layers with negative lapse rate embedded within
        inversion layers are counted as part of the inversion layer."""

        thin_layers = ((layer_df['layer_strength'] < 0) &
                       (layer_df['layer_depth'] < params['max_embed_depth'])).values
        if np.any(thin_layers):
            for x in np.argwhere(thin_layers):
                if (x > 0) & (x < len(layer_df)-1):
                    if np.all([layer_df.loc[x-1, 'layer_strength'].values > 0,
                               layer_df.loc[x+1, 'layer_strength'].values > 0]):
                        idxb = int(layer_df.loc[x, 'idx_base']) + 1
                        idxt = int(layer_df.loc[x, 'idx_top'])
                        sign_vector[idxb:idxt] = 1
        return sign_vector

    updating = True
    # If any inversions:
    layer_df = build_layer_df(sign_vector, data)
    if len(layer_df) == 1:
        updating = False
    while updating:
        sign_vector = merge_iso_layers(sign_vector, params)
        sign_vector = drop_thin_inv_layers(sign_vector, layer_df, params)
        sign_vector = merge_neg_layers(sign_vector, layer_df, params)
        
        new_layer_df = build_layer_df(sign_vector, data)

        if layer_df.shape == new_layer_df.shape:
            updating = False
        layer_df = new_layer_df
        
    layer_df['dewpoint_change'] = layer_df['dewpoint_temperature_top'] - layer_df['dewpoint_temperature_base']

    # get uncertainty of inversion strength and depth
    layer_df['uncertainty_layer_strength'] = (layer_df['uncertainty_temperature_top']**2 + 
                                              layer_df['uncertainty_temperature_base']**2)**0.5
    layer_df['uncertainty_layer_depth'] = (layer_df['uncertainty_altitude_top']**2 + 
                                           layer_df['uncertainty_altitude_base']**2)**0.5
    layer_df['uncertainty_layer_dewpoint_change'] = (layer_df['uncertainty_dewpoint_temperature_top']**2 + 
                                                     layer_df['uncertainty_dewpoint_temperature_base']**2)**0.5
    
    if params['classify_inversions']:
        layer_df['classification'] = np.nan
        layer_df.loc[layer_df['altitude_base'] == 0, 'classification'] = 'SBI'
        layer_df.loc[layer_df['altitude_base'] > 0, 'classification'] = 'EI'
        layer_df.loc[(layer_df['dewpoint_change'] > 0) & (layer_df['altitude_base'] > 0), 'classification'] = 'EI-WAA'
        layer_df.loc[(layer_df['dewpoint_change'] < 0) & (layer_df['altitude_base'] > 0), 'classification'] = 'EI-AC'    
    layer_df = layer_df[layer_df['layer_strength'] > 0].reset_index(drop=True)
    layer_df.index += 1
    return layer_df

def find_inversion_layers(dataset, params):
    """Compute the lapse rates and uncertainty. Use the k param
    to determine significance of inversion layers. k=0 corresponds 
    to methods in previous work.
    
    Note on indexing: I used forward differences, so that means a positive
    lapse rate index i is an estimate of the lapse rate from index i 
    to index i+1.
    """
    
    if 'lapse_rate' not in dataset:
        dataset = add_lapse_rate(dataset, params)
    inv_list = []
    for name, group in dsl.groupby('date'):
        inv_list.append(find_inversions(group.to_dataframe(), params))
    return pd.concat(inv_list)