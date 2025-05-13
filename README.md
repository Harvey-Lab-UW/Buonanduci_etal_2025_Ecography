# Buonanduci_etal_2025_Ecography
Data and code accompanying the manuscript 'Patterns and drivers of biotic disturbance hotspots in western United States coniferous forests' by Buonanduci, Hart, Tobin, and Harvey, published in Ecography.


[![CC BY 4.0][cc-by-shield]][cc-by]

This information is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. Any user of these data ("User" hereafter) is required to cite it appropriately in any publication that results from its use. These data may be actively used by others for ongoing research, so coordination may be necessary to prevent duplicate publication. The User is urged to contact the authors of these data for questions about methodology or results.  The User is encouraged to consider collaboration or co-authorship with authors where appropriate. Misinterpretation of data may occur if used out of context of the original study. Substantial efforts are made to ensure accuracy of the data and documentation, however complete accuracy of data sets cannot be guaranteed. All data are made available as is. Data may be updated periodically and it is the responsibility of the User to check for new versions of the data. The authors and the repository where these data were obtained shall not be liable for damages resulting from any use or misinterpretation of the data.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



For reproducibility, the following files are made available:

## Data files

#### C_input.csv
This file contains hotspots data for the Cascades region, along with covariates included in the analysis. The following columns are included:

- **damage_year**: Focal year for tree damage detection.
- **X_m**: Easting of the 5.1-km raster cell centroid (meters). Coordinate reference system is NAD 1983 USGS Contiguous USA Albers (ESRI:102039). 
- **Y_m**: Northing of the 5.1-km raster cell centroid (meters). Coordinate reference system is NAD 1983 USGS Contiguous USA Albers (ESRI:102039). 
- **Host_rich**: Host tree richness, calculated as the average host species richness within subcells containing $\geq$ 2 potential host tree species.
- **Host_BA**: Host tree basal area (m<sup>2</sup>/ha), calculated as the average total basal area of host species within subcells containing $\geq$ 2 potential host tree species.
- **Host_presence**: Prevalence or extent of host tree co-occurrence, calculated as the number of 510-m subcells containing $\geq$ 2 potential host tree species.
- **Elevation**: Elevation (m).
- **TWI**: Topographic wetness index, an index of the long-term moisture availability of a given site in the landscape.
- **HLI**: Heat load index, an index of potential direct incident radiation.
- **AET_norm**: Annual average actual evapotranspiration (AET; mm), expressed as 30-year normals (1991 – 2020).
- **tmin_norm**: Winter minimum temperature (°C). Calculated as average daily minimums for December through February, expressed as 30-year normals (1991 – 2020).
- **vpdmax_norm**: Summer maximum vapor pressure deficit (VPD; hPa). Calculated as the average daily maximums for June through August, expressed as 30-year normals (1991 – 2020).
- **tmin_anom**: Winter minimum temperature anomalies (°C). Calculated as average daily minimums for December through February, averaged over the 3-year hotspot detection window and expressed as deviations from 30-year normals.
- **vpdmax_anom**: Summer maximum VPD anomalies (hPa). Calculated as average daily maximums for June through August, averaged over the 3-year hotspot detection window and expressed as deviations from 30-year normals.
- **hotspot_bin**: Binary indicator for hotspot detection in one or more subcells (referred to in manuscript as *hotspot occurrence*). Equal to 1 if damage caused by $\geq$ 2 biotic agents and affecting $\geq$ 2 host tree species was detected within at least one 510-m subcell within a three-year window, including the focal year and previous two years.
- **hotspot_gt0**: Count of 510-m subcells in which hotspots were detected (referred to in manuscript as *hotspot prevalence*). 
- **Ntrials**: Number of subcells with the potential for hotspot detection. Subcells had the potential for hotspot detection if they were surveyed every year of the three-year hotspot window and were likely to contain $\geq$ 2 host tree species.
- **Ntrials_gt0**: Same as above, but *NA* where *hotspot_bin* = 0.

#### MR_input.csv
This file contains hotspots data for the Middle Rockies region, along with covariates included in the analysis. Columns are the same as described above for the Cascades.

#### SR_input.csv
This file contains hotspots data for the Southern Rockies region, along with covariates included in the analysis. Columns are the same as described above for the Cascades.

#### agents.csv
This file contains unique codes and common names for the 13 tree-killing Scolytinae species and three multi-agent mortality complexes included in the analysis.

#### agents_hotspots_contrib.csv
This file contains data describing the contributions of individual biotic agents and multi-species complexes to observed hotspots, quantified as proportions of observed 510-m hotspots in which each biotic agent was detected.



## Spatial data files

#### Cascades.shp
This file contains a polygon perimeter for the Cascades study region. Coordinate reference system is NAD 1983 USGS Contiguous USA Albers (ESRI:102039). 

#### M_Rockies.shp
This file contains a polygon perimeter for the Middle Rockies study region. Coordinate reference system is NAD 1983 USGS Contiguous USA Albers (ESRI:102039). 

#### S_Rockies.shp
This file contains a polygon perimeter for the Southern Rockies study region. Coordinate reference system is NAD 1983 USGS Contiguous USA Albers (ESRI:102039). 



## R script files

#### analysis_INLA.R
Code for fitting spatio-temporal models to hotspots data using R-INLA.

#### figure_study_regions.R
Code for producing Figure 1 in the manuscript.

#### figure_observed_hotspots.R
Code for producing Figures 2, 3, and 4 in the manuscript.

#### figure_agents_hotspots_contrib.R
Code for producing Figure 5 in the manuscript.

#### figure_coefficients.R
Code for producing Figure 6 in the manuscript.

#### figure_fitted_fixed_random.R
Code for producing Figure 7 in the manuscript.
