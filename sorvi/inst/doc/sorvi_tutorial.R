## ----install, eval=FALSE-------------------------------------------------
#  install.packages("sorvi")

## ----test, message=FALSE, warning=FALSE, eval=TRUE-----------------------
library(sorvi)

## ----install2, eval=FALSE------------------------------------------------
#  library(devtools)
#  install_github("ropengov/sorvi")

## ----locale, eval=FALSE--------------------------------------------------
#  Sys.setlocale(locale="UTF-8")

## ----translate, message=FALSE, eval=FALSE--------------------------------
#  translations <- load_sorvi_data("translations")
#  kable(as.matrix(translations))

## ----municipalityMML, message=FALSE, warning=FALSE, eval=FALSE-----------
#  municipality.info.mml <- get_municipality_info_mml()
#  library(knitr)
#  kable(municipality.info.mml[1:2,])

## ----province2, message=FALSE, warning=FALSE, echo=TRUE, eval=FALSE------
#  m2p <- municipality_to_province()
#  kable(head(m2p)) # Just show the first ones

## ----province6, message=FALSE, warning=FALSE, echo=TRUE, eval=FALSE------
#  municipality_to_province(c("Helsinki", "Tampere", "Turku"))

## ----province7, message=FALSE, warning=FALSE, echo=TRUE, eval=FALSE------
#  m2p <- municipality_to_province(c("Helsinki", "Tampere", "Turku"), municipality.info.mml)
#  kable(head(m2p))

## ----province3, message=FALSE, echo=TRUE, eval=FALSE---------------------
#  convert_municipality_codes(municipalities = c("Turku", "Tampere"))

## ----province4, message=FALSE, echo=TRUE, eval=FALSE---------------------
#  convert_municipality_codes(ids = c(853, 837))

## ----province5, message=FALSE, echo=TRUE, eval=FALSE---------------------
#  municipality_ids <- convert_municipality_codes()
#  kable(head(municipality_ids)) # just show the first entries

## ----sorvi-synonymes-1, message=FALSE------------------------------------
f <- system.file("extdata/municipality_synonymes.csv", package = "sorvi")
synonymes <- read.csv(f, sep = "\t")		 

## ----sorvi-synonymes-2, message=FALSE, eval=FALSE------------------------
#  synonymes <- check_synonymes(synonymes, include.lowercase = TRUE)

## ----sorvi-synonymes-3, message=FALSE, eval=FALSE------------------------
#  harmonized <- harmonize_names(c("Mantta", "Koski.Tl"), synonymes)
#  kable(harmonized)

## ----hetu, message=FALSE-------------------------------------------------
library(sorvi)
hetu("111111-111C")

## ----hetuvec, message=FALSE----------------------------------------------
library(knitr)
kable(hetu(c("010101-0101", "111111-111C")))

## ----hetuextract, message=FALSE------------------------------------------
hetu(c("010101-0101", "111111-111C"), extract = "gender")

## ----hetu2, fig.message=FALSE--------------------------------------------
valid_hetu("010101-0101") # TRUE/FALSE

## ----regressionline, message=FALSE, eval=TRUE, fig.width=10, fig.height=5----
library(sorvi) 
data(iris)
p <- regression_plot(Sepal.Length ~ Sepal.Width, iris) 
print(p)

## ----citation, message=FALSE, eval=TRUE----------------------------------
citation("sorvi")

## ----sessioninfo, message=FALSE, warning=FALSE---------------------------
sessionInfo()

