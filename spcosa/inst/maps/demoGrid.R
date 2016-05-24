library(sp)
library(maptools)
grdField <- expand.grid(s1 = 1:100, s2 = 1:50)
grdField$id <- 1
coordinates(grdField) <- ~ s1 * s2
gridded(grdField) <- TRUE
writeAsciiGrid(x = grdField, fname = "demoGrid.asc")
