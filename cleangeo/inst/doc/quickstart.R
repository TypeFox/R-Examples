## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("eblondel/cleangeo")

## ------------------------------------------------------------------------
library(cleangeo)

## ------------------------------------------------------------------------
file <- system.file("extdata", "example.shp", package = "cleangeo")

require(maptools)
sp <- readShapePoly(file)

## ------------------------------------------------------------------------
report <- clgeo_CollectionReport(sp)
clgeo_SummaryReport(report)

## ------------------------------------------------------------------------
sp.clean <- clgeo_Clean(sp)


## ------------------------------------------------------------------------
report.clean <- clgeo_CollectionReport(sp.clean)
clgeo_SummaryReport(report.clean)


## ------------------------------------------------------------------------
require(rgeos)
sapply(slot(sp.clean, "polygons"), function(x){
  gIsValid(SpatialPolygons(Srl = list(x)))
})


