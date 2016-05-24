### R code from vignette source 'jss1079.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: jss1079.Rnw:69-73
###################################################
if(!require(gstat)){install.packages("gstat"); library(gstat)}
if(!require(raster)){install.packages("raster"); library(raster)}
if(!require(plotKML)){install.packages("plotKML"); library(plotKML)}
if(!require(GSIF)){install.packages("GSIF"); library(GSIF)}


###################################################
### code chunk number 2: jss1079.Rnw:76-77
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 3: jss1079.Rnw:83-86
###################################################
rm(list=ls())
plotKML.version <- sessionInfo()[["otherPkgs"]][["plotKML"]][["Version"]]
plotKML.env(show.env = FALSE)


###################################################
### code chunk number 4: jss1079.Rnw:185-192
###################################################
library("sp")
lat = 37.423156
lon = -122.084917
name = "Google headquarters"
pnt = data.frame(name, lat, lon)
coordinates(pnt) <- ~lon+lat
proj4string(pnt) <- CRS("+proj=longlat +datum=WGS84")


###################################################
### code chunk number 5: jss1079.Rnw:213-223
###################################################
library("XML")
pnt.kml <- newXMLNode("kml")
h2 <- newXMLNode("Document", parent = pnt.kml)
h3 <- newXMLNode("name", "Google headquarters", parent = h2)
h4 <- newXMLNode("Folder", parent=pnt.kml[["Document"]])
txtc <- sprintf('<Placemark><Point><coordinates>
    %.5f,%.5f,%.0f</coordinates></Point></Placemark>',
    coordinates(pnt)[,1], coordinates(pnt)[,2], rep(0, length(pnt)))
parseXMLAndAdd(txtc, parent = h4)
pnt.kml


###################################################
### code chunk number 6: jss1079.Rnw:255-260
###################################################
library("plotKML")
data("eberg")
eberg <- eberg[runif(nrow(eberg))<.2,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")


###################################################
### code chunk number 7: jss1079.Rnw:265-266
###################################################
eberg.ll <- reproject(eberg)


###################################################
### code chunk number 8: jss1079.Rnw:271-272
###################################################
kml(eberg.ll["CLYMHT_A"], colour=CLYMHT_A)


###################################################
### code chunk number 9: jss1079.Rnw:289-291
###################################################
spplot(eberg.ll["CLYMHT_A"], edge.col="black",
  alpha=0.8, cex=seq(.3,3,length=5))


###################################################
### code chunk number 10: jss1079.Rnw:326-329
###################################################
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(eberg.ll, shape = shape, colour = CLYMHT_A, labels = SNDMHT_A,
  altitude = SNDMHT_A*10, extrude = TRUE)


###################################################
### code chunk number 11: jss1079.Rnw:366-370
###################################################
data("eberg_grid")
coordinates(eberg_grid) <- ~x+y
gridded(eberg_grid) <- TRUE
proj4string(eberg_grid) <- CRS("+init=epsg:31467")


###################################################
### code chunk number 12: jss1079.Rnw:375-379
###################################################
kml_open("eberg.kml")
kml_layer(eberg_grid, colour=TWISRT6)
kml_layer(eberg.ll[1,], colour=CLYMHT_A)
kml_close("eberg.kml")


###################################################
### code chunk number 13: jss1079.Rnw:418-420
###################################################
library("sp")
demo(meuse, echo=FALSE)


###################################################
### code chunk number 14: jss1079.Rnw:425-429
###################################################
library("GSIF")
omm <- fit.gstatModel(meuse, om~dist+ffreq, meuse.grid,
   family = gaussian(log))
om.rk <- predict(omm, meuse.grid)


###################################################
### code chunk number 15: jss1079.Rnw:475-488
###################################################
data("fmd")
fmd0 <- data.frame(fmd)
coordinates(fmd0) <- c("X", "Y")
proj4string(fmd0) <- CRS("+init=epsg:27700")
fmd_sp <- as(fmd0, "SpatialPoints")
dates <- as.Date("2001-02-18")+fmd0$ReportedDay
library("spacetime")
fmd_ST <- STIDF(fmd_sp, dates, data.frame(ReportedDay=fmd0$ReportedDay))
data("northcumbria")
ln <- Line(northcumbria)
NC <- SpatialLines(list(Lines(list(ln), ID="NC")))
proj4string(NC)  <- CRS("+init=epsg:27700")
stplot(fmd_ST, sp.layout=list("sp.lines", NC), col.regions=SAGA_pal[[1]])


###################################################
### code chunk number 16: jss1079.Rnw:495-496
###################################################
kml(fmd_ST, colour=ReportedDay, colour_scale=SAGA_pal[[1]])


