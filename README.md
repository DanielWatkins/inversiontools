# inversiontools
inversiontools is a set of functions for identifying and characterizing temperature inversions from radiosonde profile data. It was designed to be used with data from the Integrated Global Radiosonde Archive retrieved using the Siphon library, however it includes utilities for data in different formats. Temperature inversions are defined a variety of ways in meteorology and climate science. The functions here are designed to recreate common methods for quantifying temperature inversions while allowing customization (important especially for analyzing sensitivity of parameter choices). 

The process of selecting sounding layers and merging them is flexible and can be applied to different aspects of the data. The function cloud_layers merges layers based on the relative humidity.

Utilities  
- significant_levels(data, params)
- layer_average?
- classify_by_height() - classify sbi, ele, ele2, ele_sbi, etc based on position
- classify_by_type() - use the dewpoint temperatures to classify elevated inversions based on possible origin
- curvature() - estimate the curvature of the potential temperature to flag radiative inversions
- compute_uncertainty() - using either fixed or variable uncertainty estimate local lapse rate uncertainty

## Structure
The core of the code is making and updating a sign vector that marks whether or not layers are included. The merge-layers functions takes that dataframe and collapses it based on transitions in the sign vector. The process iterates:
1. Make an initial sign vector. The formation of the sign vector depends on the variable tested, thresholds, etc.
2. Form a layer dataframe.
3. Test the




## Functions
### kahl_inversions
Description of the algorithm:
> Beginning at the surface, each temperature profile is scanned upward to locate the first layer in which the temperature increases with altitude. The bottom of this layer identifies the inversion base. The inversion top is defined as the bottom of the first subsequent layer in which the temperature decreases with altitude. Thin negative-lapse layers (<100 m) are ignored if they are embedded within a deeper inversion layer.

References: Kahl, J. D. (1990). Characteristics of the low-level temperature inversion along the Alaskan Arctic coast. International Journal of Climatology, 10, 237–548.

### lower_tropospheric_stability
Lower tropospheric stability is a bulk estimate of stability used in tropical stratocumulus cloud research and in Arctic atmospheric boundary layer research. In Arctic ABL research it is often defined as the temperature difference between the 850 mb and 1000 mb pressure levels, although using 925 mb for the upper level and/or using 2-m air temperature for the lower level are not uncommon. Another variant is to use potential temperature rather than dry-bulb temperature; potential temperature differences are more important for potential for buoyant mixing.

### inversion_activity

### maximum_temperature


### Fixed geometric height
Mostly used in studies using ground stations, e.g. paired rooftop and surface station, towers, or paired valley/ridge stations.



### Inversion Activity
Defined by Milionis and Davies as the depth-weighted average difference in potential temperature across all inversion levels in a sounding.

### Maximum Temperature
Used in Connoley's analysis of surface inversions in Antarctica. First identify the height of the temperature maximum, then average those heights over a month. That results in the inversion depth, and the difference between average maximum temperature and the average surface temperature is the inversion strength.
