## ----message=FALSE-------------------------------------------------------
library(internetarchive)
library(dplyr)

## ------------------------------------------------------------------------
ia_keyword_search("isaac hecker")

## ---- eval=FALSE---------------------------------------------------------
#  ia_browse("TheLifeOfFatherHecker")

## ------------------------------------------------------------------------
ats_query <- c("publisher" = "american tract society", "year" = "1864")
ia_search(ats_query, num_results = 3)

## ---- eval=FALSE---------------------------------------------------------
#  ia_search(c("publisher" = "american tract society", date = "1840 TO 1850"))

## ---- eval=FALSE---------------------------------------------------------
#  hecker <- ia_get_items("TheLifeOfFatherHecker")

## ---- eval=FALSE---------------------------------------------------------
#  ia_metadata(hecker)
#  ia_files(hecker)

## ---- eval=FALSE---------------------------------------------------------
#  ia_keyword_search("isaac hecker", num_results = 3) %>%
#    ia_get_items() %>%
#    ia_metadata() %>%
#    filter(field == "title") %>%
#    select(value)

## ---- eval=FALSE---------------------------------------------------------
#  dir <- tempdir()
#  ia_search(ats_query, num_results = 2) %>%
#    ia_get_items() %>%
#    ia_files() %>%
#    filter(type == "txt") %>%
#    group_by(id) %>%
#    slice(1) %>%
#    ia_download(dir = dir, overwrite = FALSE) %>%
#    glimpse()

