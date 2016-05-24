library(adegraphics)
library(spdep)
library(maptools)
library(lattice)
pdf("nbgraph.pdf")

columbus <- readShapePoly(system.file("etc/shapes/columbus.shp", package = "spdep")[1])
coords <- coordinates(columbus)
col.gal.nb <- read.gal(system.file("etc/weights/columbus.gal", package = "spdep")[1])

nbobject <- col.gal.nb
xyplot(coords[, 2] ~ coords[, 1],
  		 panel = function(...) {adeg.panel.nb(col.gal.nb, coords)})

g1 <- s.label(coords, nb = nbobject, porigin.include = F, plabels.cex = 0.7, ppoints.cex = 2, Sp = columbus, pSp.col = "red", pSp.alpha = 0.5)
