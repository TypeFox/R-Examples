library(adegraphics)
library(maptools)
pdf("panelSpatial.pdf")

## ex1
nc <- readShapePoly(system.file("shapes/sids.shp", package = "maptools")[1], proj4string = CRS("+proj=longlat +datum=NAD27"))
dfxy1 <- coordinates(nc)
g1 <- s.label(dfxy1, Sp = nc, pSp.col = colorRampPalette(c("yellow", "blue"))(52), pgrid.draw = FALSE, plabels.cex = 0)


## ex2
data(meuse, package = "sp")
coordinates(meuse) <- ~ x + y
data(meuse.grid)
m <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")], data = meuse.grid)
data(meuse.riv)
meuse.sr <- SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), "meuse.riv")))

scale1 <- list("SpatialPolygonsRescale", offset = c(179900, 329600), scale = 500, fill = c("transparent", "black"), layout.scale.bar())
text11 <- list("sp.text", c(179900, 329700), "0")
text21 <- list("sp.text", c(180400, 329700), "500 m")
arrow1 <- list("SpatialPolygonsRescale", offset = c(178750, 332500), scale = 400, layout.north.arrow())
river <- list("sp.polygons", meuse.sr, fill = "lightblue")
dfxy2 <- as.data.frame(coordinates(meuse))
g2 <- s.value(dfxy2, z = meuse[, 1]$cadmium, sp.layout = list(scale1, text11, text21, arrow1, river), Sp = m)

fac <- meuse@data$ffreq
levels(fac)[1] <- "1 in 2 years"
levels(fac)[2] <- "1 in 10 years"
levels(fac)[3] <- "1 in 50 years"
arrow2 <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(181750, 330000), scale = 400)
scale2 <- list("SpatialPolygonsRescale", layout.scale.bar(), offset =  c(178050, 333600), scale = 500, fill = c("transparent", "black"))
text12 <- list("sp.text", c(178050, 333700), "0")
text22 <- list("sp.text", c(178550, 333700), "500 m")
g3 <- s.class(dfxy2, fac = fac, sp.layout = list(scale2, text12, text22, arrow2, river), starSize = 1, col = c(1, 2, 4), pellipses.col = c(1, 2, 4), 
  pellipses.alpha = 0.7, plines.lty = 3, psub.text = "Flooding frequency \n near the Meuse river", psub.pos = c(0.2, 0.88), pgrid.text.cex = 0,
  porigin.include = FALSE, Sp = meuse.sr)

## ex3
library(Guerry)
data(gfrance85)

dfxy4 <- coordinates(gfrance85)
region.names <- data.frame(gfrance85)[, 5]
col.region <- colors()[c(149, 254, 468, 552, 26)]
g4 <- s.class(dfxy4, region.names, ellip = 0, star = 0, col = col.region, Sp = gfrance85, pSp.col = col.region[region.names], porig.inclu = F)


## ex4
library(sp)
library(lattice)
nc <- readShapePoly(system.file("shapes/sids.shp", package = "maptools")[1], proj4string = CRS("+proj=longlat +datum=NAD27"))

sp <- SpatialPolygons(nc@polygons, nc@plotOrder)
g5 <- xyplot(1 ~ 1, xlim = bbox(sp)[1, ], ylim = bbox(sp)[2, ], panel = function(...) {adeg.panel.Spatial(SpObject = sp, col = "black", border = "black")})
g6 <- xyplot(1 ~ 1, xlim = bbox(sp)[1, ], ylim = bbox(sp)[2, ], panel = function(...) {adeg.panel.Spatial(nc, col = 1:14, border = "black")})
g7 <- xyplot(1 ~ 1, xlim = bbox(sp)[1, ], ylim = bbox(sp)[2, ], aspect = "iso", panel = function(...) {sp.polygons(nc, col = "black", fill = 1:5)})
g8 <- xyplot(1 ~ 1, xlim = bbox(sp)[1, ], ylim = bbox(sp)[2, ], panel = function(...) {adeg.panel.Spatial(SpObject = sp, col = "black", border = "blue")})
g9 <- xyplot(1 ~ 1, xlim = bbox(sp)[1, ], ylim = bbox(sp)[2, ], panel = function(...) {adeg.panel.Spatial(SpObject = nc, col = "black", border = "blue")})
#g10 <- s.label(cbind(-80, 35), Sp = nc)
#g11 <- s.label(cbind(-80, 35), Sp = sp)


## ex5
data(jv73, package = "ade4")
g12 <- s.label(jv73$xy, Sp = jv73$Spatial)
g13 <- s.label(jv73$xy, Sp = jv73$Spatial, pSp.col = "red")

spoints <- SpatialPoints(jv73$xy)
g14 <- s.label(jv73$xy, Sp = spoints, plab.cex = 0, ppoin.cex = 0, pSp.col = 1)

sgrid <- SpatialGrid(GridTopology(c(0, 0), c(1, 1), c(3, 5)))
xyplot(0:5 ~ 0:3, panel = function(...) sp.grid(sgrid, col = 1:2))

nc <- SpatialGridDataFrame(getGridTopology(sgrid), data = data.frame(matrix(1:15, ncol = 1)))
xyplot(0:5 ~ 0:3, panel = function(...) sp.grid(nc, col = 1, at = pretty(rnorm(15), 2), col.region = 2:3))
xyplot(0:5 ~ 0:3, panel = function(...) adeg.panel.Spatial(nc, col = 1:3))
xyplot(0:5 ~ 0:3, panel = function(...) adeg.panel.Spatial(nc, col = 1:2))


## ex6
mysp <- SpatialPointsDataFrame(matrix(rnorm(20), 10), data.frame(matrix(rnorm(20), 10)))
s.Spatial(mysp, pSp.cex = 2)
s.Spatial(mysp, col = c("red", "blue"), pSp.cex = 2)
s.Spatial(mysp, pSp.col = c("red", "blue"), pSp.cex = 2) ## must be the same as the previous plot
