## ----install1, eval=FALSE------------------------------------------------
#  install.packages("pxweb")

## ----install2, eval=FALSE------------------------------------------------
#  library("devtools")
#  devtools::install_github("ropengov/pxweb")

## ----test, message=FALSE, warning=FALSE, eval=TRUE-----------------------
library(pxweb)

## ----locale, eval=FALSE--------------------------------------------------
#  Sys.setlocale(locale="UTF-8")

## ----standardquery, message=FALSE, eval=FALSE----------------------------
#  # Navigate through all pxweb api:s installed.
#  d <- interactive_pxweb()
#  
#  # Get data from SCB (Statistics Sweden)
#  d <- interactive_pxweb(api = "api.scb.se")
#  
#  # Fetching data from the swedish SCB (Statistics Sweden) pxweb API:
#  d <- interactive_pxweb(api = "api.scb.se", version = "v1", lang = "sv")
#  
#  # Fetching data from statfi (Statistics Finland)
#  d <- interactive_pxweb(api = "pxwebapi2.stat.fi")

## ----directquery, message=FALSE, eval=FALSE------------------------------
#  pxweb_test_data <-
#    get_pxweb_data(url = "http://api.scb.se/OV0104/v1/doris/sv/ssd/PR/PR0101/PR0101E/Basbeloppet",
#                   dims = list(ContentsCode = c('PR0101A1'),
#                               Tid = c('*')),
#                   clean = FALSE)

## ----citation, message=FALSE, eval=TRUE----------------------------------
citation("pxweb")

## ----sessioninfo, message=FALSE, warning=FALSE---------------------------
sessionInfo()

