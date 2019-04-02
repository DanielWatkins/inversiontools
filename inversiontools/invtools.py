# Copyright (c) 2018 Daniel Watkins

"""Functions to compute inversions following the methodologies used in Arctic ABL research.
"""


def inv_finder(data,min_lapse,min_depth,iso_base, iso_adj, max_embed_depth, return_mode):
    """Locate the base and top of each inversion in a profile.

    Inputs
    ------
    data : Pandas DataFrame 
    mode : string
        Options: "scanning", "pressure", "max_t".
        scanning : use the full profile.
        pressure : use only the difference in temp. between two specified pressures
        max_t : return the difference between maximum temp. in a profile and the surface temp.

    min_lapse : float
        minimum dt/dz (deg C per km) to consider a layer to be an inversion layer.
        Default 0 means any positive lapse rate counts as an inversion layer.

    max_lapse : float
        maximum dt/dz (deg C/km) for a sounding to be considered valid.
        Default is 500, which is equivalent to 50 degrees in 100 meters. That's
        fairly high, but since most studies don't use this parameter the default
        is set high enough to let most soundings through.

    min_depth : float
        minimum depth of an inversion layer (m). Layers with depth less than
        min_depth will be discarded. Default is 20 m, based on [1].

    max_embed_depth : float
        embedded noninversion layers with thickness less than max_embed_depth (m)
        will be merged. Default 100 meters, which is standard in the literature [2].

    Handling isothermal layers:

    iso_base : bool
        Consider an isothermal layer adjacent and below an inversion layer
        to be part of that layer. If iso_all = True and iso_base = False,
        inversions are required to start with a strict inversion.

    iso_adj : bool
        Consider isothermal layers adjacent to strict inversions (embedded or above)
        to be part of an inversion. Iso_base is considered separately.
    
    pressure_levels : int tuple
        Pressure levels to use for differencing. Default is (850,1000), e.g. [5].

    return_mode:     Options: "all","sbi","strongest","first"                
                     "all" means return the dataframe with all inversions after merging
                     "sbi" means returning only (strictly) surface-based inversions
                     "main" means returning the inversion layer with greatest dt/dz.                     
    Returns:
    ------
    Pandas DataFrame with detected inversions.
    """

    import numpy as np
    import pandas as pd
    
    def merge_layers(idx1, idx2, layer_array):
        """Take the layers layer_array[idx1,:] and
        layer_array[idx2,:], and create a new array with
        layer_array[idx1:idx2, :] as a single layer."""
        nr, nc = layer_array.shape
        new_array = np.zeros((nr-idx2+idx1,nc))
        if idx1 > 0:
            new_array[0:idx1,:] = layer_array[0:idx1,:]

        new_array[idx1,0] = layer_array[idx1,0]
        new_array[idx1,1] = layer_array[idx2,1]
        new_array[idx1,2] = layer_array[idx1,2]
        new_array[idx1,3] = layer_array[idx2,3]

        if idx1 < nr-idx2+idx1:
            new_array[idx1+1:,:] = layer_array[idx2+1:,:]
    
        new_array[:,4] = new_array[:,1] - new_array[:,0]
        new_array[:,5] = new_array[:,3] - new_array[:,2]
        return new_array
    
    
    
    dt = np.diff(data.t.values)
    dz = np.diff(data.gph.values)
    date = data.date.values[0]
    
    # Convert min_lapse to degrees per meter not per km
    min_lapse = min_lapse/1000 

    # Lapse rates
    dtdz = dt/dz
    
    # Default inversion matrix
    df = pd.DataFrame({'date':date,'dT':np.nan,'dZ':np.nan,'bh':np.nan,'inv_num':np.nan},index=[0])
    
    
    # If no layer has a lapse rate greater than the specified minimum,
    # set inv number to 0 to show that the observation is valid, but no inv
    # is present.
    if np.all(dtdz < min_lapse):
        df['inv_num'] = 0
        return df
  
    # A sign vector will be used to determine whether inversion is present.
    sgn = np.zeros(len(dtdz))
    
    # Flag inversion / noninversion. 
    # Layers exactly equal to min lapse are considered noninversion.
    if np.any(dtdz <= min_lapse):
        sgn[dtdz <= min_lapse] = -1
    if np.any(dtdz > min_lapse):
        sgn[dtdz > min_lapse] = 1
    if np.any(dtdz == 0):
        sgn[dtdz == 0] = 0
        
    # Appending -1 to the beginning allows the inversion finder to identify
    # surface based inversions
    sgn = np.concatenate(([-1],sgn))
    
    
    # Create a list of all changes to the sign. These are the tops and bottoms of layers.
    idx=0
    idx_list=[]
    while idx < len(sgn)-1:
        if sgn[idx+1] != sgn[idx]:
            idx_list.append(idx)
        idx += 1
    
    if np.mod(len(idx_list),2)==1:
        idx_list.append(int(len(data)-1))

    # print(idx_list)
    # Store info in array. Columns:
    # Zb, Zt, Tb, Tt, dZ, dT
    nlayers = int(len(idx_list)-1)
    layer_array = np.zeros((nlayers,6))
    for i in range(0,nlayers):
        idx_b = idx_list[i]
        idx_t = idx_list[i+1]
        layer_array[i,0] = data.loc[idx_b,'gph']
        layer_array[i,1] = data.loc[idx_t,'gph']
        layer_array[i,2] = data.loc[idx_b,'t']
        layer_array[i,3] = data.loc[idx_t,'t']
        
    layer_array[:,4] = layer_array[:,1]-layer_array[:,0]
    layer_array[:,5] = layer_array[:,3]-layer_array[:,2]
    
    
    # First remove super thin layers.
    # This is to comply with the description in T&G2009
    if np.any(layer_array[:,4] < min_depth):
        idx=0
        n = layer_array.shape[0]
        while idx < n:
            if (layer_array[idx,4] < min_depth) & (layer_array[idx,4] > 0):
                layer_array = np.delete(layer_array,(idx),axis=0)
                n = layer_array.shape[0]
                idx = 0
            else:
                idx += 1
    
    if layer_array.shape[0]==0:
        return df
    
    # Join adjacent isothermal layers and remove isolated isothermal layers.
    if np.any(layer_array[:,5]==0):
        # Remove isolated isothermal layers. This section can be made to allow 
        # isothermal layers by placing it within an if statement and adding an
        # argument "iso_lone" to the function call.
        idx=0
        n = layer_array.shape[0]
        while idx < n-1:
            if (layer_array[idx,5]==0):
                if idx > 0:
                    if (layer_array[idx-1,5] < 0) & (layer_array[idx+1,5] < 0):
                        layer_array = np.delete(layer_array,(idx),axis=0)
                        idx=0
                        n = layer_array.shape[0]
                    else:
                        idx += 1
                else:
                    if layer_array[idx+1,5] < 0:
                        layer_array = np.delete(layer_array,(idx),axis=0)
                        idx=0
                        n = layer_array.shape[0]
                    else:
                        idx += 1


            else:
                idx += 1

        # If iso_base is True, then merge base isothermal layers with overlying inversion layers
        if iso_base:
            idx=0
            n = layer_array.shape[0]
            while idx < n-1:
                if (layer_array[idx,5]==0) & (layer_array[idx+1,5] > 0):
                    layer_array = merge_layers(idx, idx+1, layer_array)
                    n = layer_array.shape[0]
                    idx=0
                else:
                    idx += 1
        # If iso_adj is true, merge layers that overlie or are embedded within inversion layers.
        if iso_adj:
            idx=0
            n = layer_array.shape[0]
            # Merge sandwiched isothermal layers 
            while idx < n-2:
                if ((layer_array[idx,5]>0) & (layer_array[idx+1,5] == 0)) & (layer_array[idx+2,5]>0):
                    layer_array = merge_layers(idx, idx+2, layer_array)
                    n = layer_array.shape[0]
                    idx=0
                else:
                    idx += 1
           
            # Merge top isothermal layers. The reason this is done separately from the sandwich
            # layers is that otherwise you would get two inversion layers next to each other that should have
            # been merged but weren't.
            idx=0
            n = layer_array.shape[0]
            while idx < n-1:
                if (layer_array[idx,5]>0) & (layer_array[idx+1,5] == 0):
                    layer_array = merge_layers(idx, idx+1, layer_array)
                    n = layer_array.shape[0]
                    idx=0
                else:
                    idx += 1
    
    # Merge thin embedded non-inversion layers, delete thick embedded layers
    if np.any((layer_array[:,5] <= 0) & (layer_array[:,4] < max_embed_depth)):
        idx = 0
        n = layer_array.shape[0]
        while idx < n-2:
            if (layer_array[idx,5]>0) & (layer_array[idx+2,5]>0):
                if (layer_array[idx+1,4] < max_embed_depth) & (layer_array[idx+1,5] < 0):
                    layer_array = merge_layers(idx, idx+2, layer_array)
                    n = layer_array.shape[0]
                    idx = 0
                elif layer_array[idx+1,5] <= 0:
                    layer_array = np.delete(layer_array,(idx+1),axis=0)
                    n = layer_array.shape[0]
                    idx = 0
                else:
                    idx += 1
            else:
                idx += 1

    # Delete any remaining negative layers
    if np.any(layer_array[:,5] < 0):
        idx = 0
        n = layer_array.shape[0]
        while idx < n:
            if layer_array[idx,5] < 0:
                layer_array = np.delete(layer_array,(idx),axis=0)
                n = layer_array.shape[0]
                idx = 0
            else:
                idx += 1

    n = layer_array.shape[0]
    if n==0:
        return df
    # Finally convert into a Pandas data frame
    
    
    
    # Return according to return mode
    if return_mode == 'all':
        df = pd.DataFrame({'date':[date]*n,
                       'dT':layer_array[:,5],
                       'dZ':layer_array[:,4],
                       'bh':layer_array[:,0],
                       'inv_num':np.arange(n)+1},
                      index=np.arange(n))
        return df
    elif return_mode == 'first':
        df['dT'] = layer_array[0,5]
        df['dZ'] = layer_array[0,4]
        df['bh'] = layer_array[0,0]
        df['inv_num'] = 1
        return df
    elif return_mode == 'sbi':
        if layer_array[0,0] == data.loc[0,'gph']:
            df['dT'] = layer_array[0,5]
            df['dZ'] = layer_array[0,4]
            df['bh'] = layer_array[0,0]
            df['inv_num'] = 1
            return df
        else:
            df['inv_num'] = 0
            return df
    elif return_mode == 'main':
        if len(layer_array)==0:
            df['inv_num']=0
            return df
        elif layer_array.shape[0]==1:
            idx=0
        else:
            # Find the layer containing the maximum temperature gradient
            # Possible that this layer was deleted - what then?
            z_max_grad = data.gph[0:-1][dtdz == np.max(dtdz)].values
           
            idx = 0
            while (idx < layer_array.shape[0]) and (layer_array[idx+1,0] < z_max_grad):
                idx += 1

            
        df['dT'] = layer_array[idx,5]
        df['dZ'] = layer_array[idx,4]
        df['bh'] = layer_array[idx,0]
        df['inv_num'] = 1
        return df
        
        
    else:
        print("Bad return mode, returning all")
        df = pd.DataFrame({'date':[date]*n,
                       'dT':layer_array[:,5],
                       'dZ':layer_array[:,4],
                       'bh':layer_array[:,0],
                       'inv_num':np.arange(n)+1},
                      index=np.arange(n))
        return df

# What if the only inversion present is thinner than 20 meters?
