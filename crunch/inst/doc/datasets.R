## ---- results='hide', echo=FALSE, message=FALSE-----------------------------------------------------------------------
## Because the vignette tasks require communicating with a remote host,
## we do all the work ahead of time and save a workspace, which we load here.
## We'll then reference saved objects in that as if we had just retrieved them
## from the server
library(crunch)
load("vignettes.RData")
options(width=120)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  load("../vignettes/economist.RData")
#  dim(economist)

## ---- echo=FALSE------------------------------------------------------------------------------------------------------
dim.ds

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  ds <- newDataset(economist, name="Economist/YouGov Weekly Survey")
#  dim(ds)

## ---- echo=FALSE------------------------------------------------------------------------------------------------------
dim.ds

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  listDatasets()

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  ds <- loadDataset("Economist/YouGov Weekly Survey")

## ---------------------------------------------------------------------------------------------------------------------
is.dataset(ds)

## ---------------------------------------------------------------------------------------------------------------------
name(ds)
description(ds)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  description(ds) <- "U.S. nationally representative sample, 1000 respondents"
#  description(ds)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  ds <- refresh(ds)
#  description(ds)

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  ## Not run
#  delete(ds)

