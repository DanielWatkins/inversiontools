# inversiontools
inversiontools is a set of functions for identifying and characterizing temperature inversions from radiosonde profile data. Temperature inversions are defined a variety of ways in meteorology and climate science. The functions here are designed to recreate common methods for quantifying temperature inversions while allowing customization (important especially for analyzing sensitivity of parameter choices).

Although the functions will work on any Pandas dataframe with columns for temperature, pressure, and height, I have designed the code with the Siphon library in mind, which has tools for importing radiosonde profiles from the IGRA2 archive and from the University of Wyoming archive.

## Functions
### kahl1990
Description of the algorithm:
> Beginning at the surface, each temperature profile is scanned upward to locate the first layer in which the temperature increases with altitude. The bottom of this layer identifies the inversion base. The inversion top is defined as the bottom of the first subsequent layer in which the temperature decreases with altitude. Thin negative-lapse layers (<100 m) are ignored if they are embedded within a deeper inversion layer.



Reference: Kahl, J. D. (1990). Characteristics of the low-level temperature inversion along the Alaskan Arctic coast. International Journal of Climatology, 10, 237–548.

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
