## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)




## ----message=FALSE-------------------------------------------------------
library(rgdal)
library(raster)
library(rworldmap)
data(countriesLow)
llproj <- projection(countriesLow)
library(graticule)
map<- subset(countriesLow, SOVEREIGNT == "Australia")

## VicGrid
prj <- "+proj=lcc +lat_1=-36 +lat_2=-38 +lat_0=-37 +lon_0=145 +x_0=2500000 +y_0=2500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

pmap <- spTransform(map, CRS(prj))

## specify exactly where we want meridians and parallels
lons <- seq(140, 150, length = 5)
lats <- seq(-40, -35, length = 6)
## optionally, specify the extents of the meridians and parallels
## here we push them out a little on each side
xl <-  range(lons) + c(-0.4, 0.4)
yl <- range(lats) + c(-0.4, 0.4)
## build the lines with our precise locations and ranges
grat <- graticule(lons, lats, proj = prj, xlim = xl, ylim = yl)
## build the labels, here they sit exactly on the western and northern extent
## of our line ranges
labs <- graticule_labels(lons, lats, xline = min(xl), yline = max(yl), proj = prj)

## set up a map extent and plot
op <- par(mar = rep(0, 4))
plot(extent(grat) + c(4, 2) * 1e5, asp = 1, type = "n", axes = FALSE, xlab = "", ylab = "")
plot(pmap, add = TRUE)
## the lines are a SpatialLinesDataFrame
plot(grat, add = TRUE, lty = 5, col = rgb(0, 0, 0, 0.8))
## the labels are a SpatialPointsDataFrame, and islon tells us which kind
text(subset(labs, labs$islon), lab = parse(text = labs$lab[labs$islon]), pos = 3)
text(subset(labs, !labs$islon), lab = parse(text = labs$lab[!labs$islon]), pos = 2)
par(op)


## ----message=FALSE-------------------------------------------------------
library(raster)
library(graticule)
library(rgdal)
tfile <- system.file("extdata",  "nt_20140320_f17_v01_s.bin", package = "graticule")
ice <- raster(tfile)

meridians <- seq(-180, 160, by = 20)
parallels <- c(-80, -73.77, -68, -55, -45)
mlim <- c(-180, 180)
plim <- c(-88, -50)
grat <- graticule(lons = meridians, lats = parallels, xlim = mlim, ylim = plim, proj = projection(ice))
labs <- graticule_labels(meridians, parallels, xline = -45, yline = -60, proj = projection(ice))
plot(ice, axes = FALSE)
plot(grat, add = TRUE, lty = 3)
text(labs, lab = parse(text= labs$lab), col= c("firebrick", "darkblue")[labs$islon + 1], cex = 0.85)
title(sprintf("Sea ice concentration %s", gsub(".bin", "", basename(tfile))), cex.main = 0.8)
title(sub = projection(ice), cex.sub = 0.6)

## ------------------------------------------------------------------------
polargrid <- graticule(lons = c(meridians, 180), lats = parallels,  proj = projection(ice), tiles = TRUE)
centroids <- project(coordinates(polargrid), projection(ice), inv = TRUE)
labs <- graticule_labels(meridians, parallels,  proj = projection(ice))
labs <- graticule_labels(as.integer(centroids[,1]), as.integer(centroids[,2]),  proj = projection(ice))
labs <- labs[!duplicated(as.data.frame(labs)), ] ## this needs a fix
cols <- sample(colors(), nrow(polargrid))
op <- par(mar = rep(0, 4))
plot(polargrid, col  = cols, bg = "black")
text(labs[labs$islon, ], lab = parse(text = labs$lab[labs$islon]), col = "black",  cex = .55, pos = 3)
text(labs[!labs$islon, ], lab = parse(text = labs$lab[!labs$islon]), col = "black", cex = .55, pos = 1)
par(op)


## ------------------------------------------------------------------------
plot(pmap)
llgridlines(pmap)

## ------------------------------------------------------------------------
plot(pmap)
lons <- seq(140, 150, length = 5)
lats <- seq(-40, -35, length = 6)
llgridlines(pmap, easts = lons, norths = lats)

## ------------------------------------------------------------------------
op <- par(xpd = NA)
ex <- as(extent(range(lons), range(lats)), "SpatialPolygons")
projection(ex) <- llproj
pex <- spTransform(ex, CRS(projection(pmap)))
plot(extent(pex), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
plot(pmap, add = TRUE)
llgridlines(pex, easts = lons, norths = lats)
par(op)


## ------------------------------------------------------------------------
plot(ice, axes = FALSE)
##llgridlines(ice)  does not understand a raster
llgridlines(as(ice, "SpatialPoints"))

## ------------------------------------------------------------------------
plot(ice, axes = FALSE)
llgridlines(as(ice, "SpatialPoints"), easts = c(-180, -120, -60, 0, 60, 120), norths = c(-80, -70, -60, -50), ndiscr = 50)

## ----message=FALSE-------------------------------------------------------
library(oce)
## we need to hop the crevasse into another world
pp <- coordinates(as(as(map, "SpatialLines"), "SpatialPoints"))
mapPlot(pp[,1], pp[,2], projection = projection(pmap), longitudelim = xl, latitudelim = yl, type = "n", grid = FALSE)
mapGrid(longitude = lons, latitude = lats)
## and to prove that all is well in the world
plot(pmap, add = TRUE)


## ----message = FALSE-----------------------------------------------------
ipts <- coordinates(spTransform(xyFromCell(ice, sample(ncell(ice), 1000), spatial = TRUE), CRS(llproj)))
mapPlot(ipts[,1], ipts[,2], projection = projection(ice), type = "n", grid = FALSE)

plot(ice, add = TRUE)
mapGrid(10, 15)
mapTissot()


## ------------------------------------------------------------------------
devtools::session_info()

