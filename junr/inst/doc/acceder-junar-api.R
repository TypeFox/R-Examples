## ------------------------------------------------------------------------
library(junr)
url_base <- "http://api.datosabiertos.presidencia.go.cr/api/v2/datastreams/"
api_key <- "0bd55e858409eefabc629b28b2e7916361ef20ff" 

## ---- eval=FALSE---------------------------------------------------------
#  get_index(url_base, api_key)

## ------------------------------------------------------------------------
list_guid(url_base, api_key)

## ------------------------------------------------------------------------
list_titles(url_base, api_key)

## ------------------------------------------------------------------------
guid_datos <- "COMPR-PUBLI-DEL-MINIS"
datos_compras <- get_data(url_base, api_key, guid_datos)

## ------------------------------------------------------------------------
get_dimensions(url_base, api_key)

## ------------------------------------------------------------------------
datos_con_divisas <- get_data(url_base, api_key, "LICIT-ADJUD-POR-LOS-MINIS")

