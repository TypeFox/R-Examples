## ----setup, include=FALSE------------------------------------------------

suppressWarnings({
  suppressPackageStartupMessages({
    loadNamespace("knitr") # for opts_chunk only
    library(icd)
    library(magrittr)
    })
  })

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)

## ----vermont-charlson----------------------------------------------------
# typical hospital format data, with many columns for diagnoses
head(vermont_dx)
# convert to long format (could use other tools, but the icd version accounts better for known structure of the data.
head(vermont_dx %>% icd_wide_to_long)
# calculate charlson scores and summarize
vermont_dx %>% icd_wide_to_long %>% icd_charlson %>% summary
# show the first few actual scores: the names are the patient IDs:
vermont_dx %>% icd_wide_to_long %>% icd_charlson %>% head(25) -> vermont_charlson
vermont_charlson
names(vermont_charlson)
unname(vermont_charlson)

