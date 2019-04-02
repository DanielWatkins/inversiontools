# inversiontools
This package contains tools for identifying temperature inversions within radiosonde profiles. Temperature inversions are defined as regions where temperature increases with height. Algorithms for identifying inversions can broadly be devided into fixed level differencing methods, profile methods, and other. 

## Fixed level
### Lower Tropospheric Stability
Lower tropospheric stability is a bulk estimate of stability used in tropical stratocumulus cloud research and in Arctic atmospheric boundary layer research. It is generally defined as 
$$ LTS = T_{850 mb} - T_{1000 mb} $$
although the specific heights vary. Other options 925 mb for the upper level and SST or 2-m air temperature for the lower level.

### Fixed geometric height
Mostly used in studies using ground stations, e.g. paired rooftop and surface station, towers, or paired valley/ridge stations.

## Profile search
### Kahl's algorithm
Starting from the surface, move upward until the temperature increases. Mark that point as the inversion base. Continue until temperature decreases with height for more than 100 m. Mark that point as the inversion top. Return inversion depth (height at top minus height at base) and inversion strength (temperature at top minus temperature at base). Zhang and Seidel tested the sensitivity of this algorithm to differences in the threshold for layers of decreasing temperature and found the algorithm to not be very sensitive. 

### Inversion Activity
Defined by Milionis and Davies as the depth-weighted average difference in potential temperature across all inversion levels in a sounding.

### Maximum Temperature
Used in Connoley's analysis of surface inversions in Antarctica. First identify the height of the temperature maximum, then average those heights over a month. That results in the inversion depth, and the difference between average maximum temperature and the average surface temperature is the inversion strength.