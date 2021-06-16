# GIS Data
GIS data was assembled based on data in the CSV files, "FOCB Monitoring Sites
SHORT NAMES.xlsx" and "din_2019.csv" and "tn_18_19.csv". Those files were
generated via the R Notebooks, "FOCB_DIN_Analysis.Rmd" and
"FOCB_TN_Analysis.Rmd". See the Notebooks for details.

## FOCB Data
### Locations
Data from "FOCB Monitoring Sites SHORT NAMES.xlsx" was loaded into GIS as a
table.  A feature layer was created by displaying the XY data (longitudes and
latitudes) from that file, and exporting as a feature layer.  The data layer was
saved in a Personal Geodatabase because shapefiles generally do not accept NULL
values in their attributes.

### Attributes
Specific data on DIN and TN were loaded from "din_2019.csv" and "tn_18_19.csv".
These file contains descriptive statistics by site derived from specific subsets
of recent FOCB data. Analysis shows that year to year and seasonal variation is
substantial, so we focused on periods of time when FOCB sampling was reasonably
consistent across sites.  See the R Notebooks "FOCB_DIN_Analysis.Rmd" and 
"FOCB_TN_Analysis.Rmd" for details.

### Combining the Data
Data from the attribute tables was added to the features by importing the 
attribute data as separate table, joining the tables by Site ID,  copying 
data from the joined tables into the feature data, and finally removing the
joins.  Data imported to the feature layer represent median values for
purposes of GIS mapping.

### Attribute Definitions
Data in the dep_locations feature data have the following attributes:

| Site        |   Definition                                            |
|:------------|:--------------------------------------------------------|
|Station_ID   |   Site Code from DEP data, with depth designations (if any) removed |
|Station_Name |   Original name of the site from the DEP source data    |
|Alt_Name     |   Simplified name for graphics and GIS display          |
|Town         |   Town where the site is located
|Y            |   Latitude.  Interpreted here as WGS 1984               |
|X            |   Longitude.  Interpreted here as WGS 1984              |
|Category     |   Type of monitoring location -- Surface, Profile, or Surface/Continuous  |  
|DIN          |   Dissolved Inorganic Nitrogen (sum of NH4 and NOX values) from 2019, mg/l, median values |
|TN           |   Total nitrogen,  from 2018 and 2018 sampling, mg/l, median values |

## DEP Data
The DEP data presented here was developed in the separate analysis of DEP data,
in the "DEP_nutrient" repository at https://github.com/ccb60/DEP_nutrients.
Procedures are similar to what we did with the FOCB dat.  See that repository
for details.
