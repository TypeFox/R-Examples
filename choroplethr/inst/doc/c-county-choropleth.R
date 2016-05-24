## ----hold=TRUE-----------------------------------------------------------
library(choroplethr)

?df_pop_county
data(df_pop_county)

?county_choropleth
county_choropleth(df_pop_county)

## ------------------------------------------------------------------------
library(choroplethrMaps)

?county.regions
data(county.regions)
head(county.regions)

## ------------------------------------------------------------------------
county_choropleth(df_pop_county,
                 title      = "2012 Population Estimates",
                 legend     = "Population",
                 num_colors = 1,
                 state_zoom = c("california", "washington", "oregon"))

## ------------------------------------------------------------------------
# FIPS codes for Alameda, Contra Costa, Marin, Napa, San Francisco, San Mateo, Santa Clara, 
# Solano, and Sonoma counties
bay_area_counties = c(6001, 6013, 6041, 6055, 6075, 6081, 6085, 6095, 6097)
county_choropleth(df_pop_county,
                 title       = "2012 Population Estimates",
                 legend      = "Population",
                 num_colors  = 1,
                 county_zoom = bay_area_counties)

## ------------------------------------------------------------------------
library(ggplot2)

choro = CountyChoropleth$new(df_pop_county)
choro$title = "2012 Population Estimates"
choro$ggplot_scale = scale_fill_brewer(name="Population", palette=2, drop=FALSE)
choro$render()

