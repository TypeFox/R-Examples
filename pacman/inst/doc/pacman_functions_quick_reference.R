## ----setup, include=FALSE------------------------------------------------
# set global chunk options
library(knitr); library(pacman); library(methods)
opts_chunk$set(cache=FALSE, comment=NA, warning=FALSE)

## Function for embedding high qual text images:
uri_embed <- function(path, add="") {
    uri <- knitr::image_uri(path)
    cat(paste0("<img ", add, " src=\"", uri, "\" />"))
}
opts_knit$set(upload.fun = image_uri, self.contained=TRUE)

## set mirror
options(repos="http://cran.rstudio.com/")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
uri_embed("r_pacman.png", 
    "width=\"350\", height=\"150\" style=\"display:block; margin-left:auto; margin-right:auto;\"")

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/01_installing_loading_deleting.R")
cat(paste(installing_tab, collapse="\n"))

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/02_session_information.R")
cat(paste(installing_tab, collapse="\n"))

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/03_local_package_information.R")
cat(paste(installing_tab, collapse="\n"))

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/04_internet_based_info.R")
cat(paste(installing_tab, collapse="\n"))

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/05_library_information.R")
cat(paste(installing_tab, collapse="\n"))

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
installing_tab <- readLines("tables/06_pacman_tools.R")
cat(paste(installing_tab, collapse="\n"))

## ----clean-up, include=FALSE---------------------------------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))

