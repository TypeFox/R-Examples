## ------------------------------------------------------------------------
library(junr)
base_url <- "http://api.datosabiertos.presidencia.go.cr/api/v2/datastreams/"
api_key <- "0bd55e858409eefabc629b28b2e7916361ef20ff" 

## ---- eval=FALSE---------------------------------------------------------
#  get_index(base_url, api_key)

## ------------------------------------------------------------------------
list_guid(base_url, api_key)

## ------------------------------------------------------------------------
pres_list <-list_guid(base_url, api_key)
pres_list[3]

## ------------------------------------------------------------------------
list_titles(base_url, api_key)

## ------------------------------------------------------------------------
data_guid <- "COMPR-PUBLI-DEL-MINIS"
purchasing_data <- get_data(base_url, api_key, data_guid)

## ------------------------------------------------------------------------
get_dimensions(base_url, api_key)

## ------------------------------------------------------------------------
currency_data <- get_data(base_url, api_key, "LICIT-ADJUD-POR-LOS-MINIS")
currency_data$`Monto Adjudicado` <- clean_currency(currency_data$`Monto Adjudicado`)  

