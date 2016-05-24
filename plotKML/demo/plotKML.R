## Complete tutorial available at: [http://gsif.isric.org/doku.php?id=wiki:tutorial_plotkml]

plotKML.env(silent = FALSE, kmz = FALSE)
## -------------- SpatialPointsDataFrame --------- ##
library(sp)
library(rgdal)
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
## subset to 20 percent:
eberg <- eberg[runif(nrow(eberg))<.1,]
## bubble type plot:
plotKML(eberg["CLYMHT_A"])
plotKML(eberg["CLYMHT_A"], colour_scale=rep("#FFFF00", 2), points_names="")

## -------------- SpatialLinesDataFrame --------- ##
data(eberg_contours)
plotKML(eberg_contours)
## plot contour lines with actual altitudes:
plotKML(eberg_contours, colour=Z, altitude=Z)

## -------------- SpatialPolygonsDataFrame --------- ##
data(eberg_zones)
plotKML(eberg_zones["ZONES"])
## add altitude:
zmin = 230
plotKML(eberg_zones["ZONES"], altitude=zmin+runif(length(eberg_zones))*500)

## -------------- SpatialPixelsDataFrame --------- ##
library(rgdal)
library(raster)
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
TWI <- reproject(eberg_grid["TWISRT6"])
data(SAGA_pal)
plotKML(TWI, colour_scale = SAGA_pal[[1]])
## set limits manually (increase resolution):
plotKML(TWI, z.lim=c(12,20), colour_scale = SAGA_pal[[1]],
  png.width = gridparameters(TWI)[1,"cells.dim"]*5, 
  png.height = gridparameters(TWI)[2,"cells.dim"]*5)

## categorical data:
eberg_grid$LNCCOR6 <- as.factor(paste(eberg_grid$LNCCOR6))
levels(eberg_grid$LNCCOR6)
data(worldgrids_pal)
## attr(worldgrids_pal["corine2k"][[1]], "names")
pal = as.character(worldgrids_pal["corine2k"][[1]][c(1,11,13,14,16,17,18)])
LNCCOR6 <- reproject(eberg_grid["LNCCOR6"])
plotKML(LNCCOR6, colour_scale=pal)

## -------------- SpatialPhotoOverlay --------- ##
library(RCurl)
imagename = "Soil_monolith.jpg"
urlExists = url.exists("http://commons.wikimedia.org")
if(urlExists){
  x1 <- getWikiMedia.ImageInfo(imagename)
  sm <- spPhoto(filename = x1$url$url, exif.info = x1$metadata)
  # str(sm)
  plotKML(sm)
}

## -------------- SoilProfileCollection --------- ##
library(aqp)
library(plyr)
## sample profile from Nigeria:
lon = 3.90; lat = 7.50; id = "ISRIC:NG0017"; FAO1988 = "LXp" 
top = c(0, 18, 36, 65, 87, 127) 
bottom = c(18, 36, 65, 87, 127, 181)
ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
hue = c("7.5YR", "7.5YR", "2.5YR", "5YR", "5YR", "10YR")
value = c(3, 4, 5, 5, 5, 7); chroma = c(2, 4, 6, 8, 4, 3)
## prepare a SoilProfileCollection:
prof1 <- join(data.frame(id, top, bottom, ORCDRC, hue, value, chroma), 
   data.frame(id, lon, lat, FAO1988), type='inner')
prof1$soil_color <- with(prof1, munsell2rgb(hue, value, chroma))
depths(prof1) <- id ~ top + bottom
site(prof1) <- ~ lon + lat + FAO1988 
coordinates(prof1) <- ~ lon + lat
proj4string(prof1) <- CRS("+proj=longlat +datum=WGS84")
prof1
plotKML(prof1, var.name="ORCDRC", color.name="soil_color")

## -------------- STIDF --------- ##
library(spacetime)
## daily temperatures for Croatia:
data(HRtemp08)
## format the time column:
HRtemp08$ctime <- as.POSIXct(HRtemp08$DATE, format="%Y-%m-%dT%H:%M:%SZ")
## create a STIDF object:
sp <- SpatialPoints(HRtemp08[,c("Lon","Lat")])
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84")
HRtemp08.st <- STIDF(sp, time = HRtemp08$ctime, data = HRtemp08[,c("NAME","TEMP")])
## subset to first 500 records:
HRtemp08_jan <- HRtemp08.st[1:500]
str(HRtemp08_jan)
plotKML(HRtemp08_jan[,,"TEMP"], dtime = 24*3600, LabelScale = .4)

## foot-and-mouth disease data:
data(fmd)
fmd0  <- data.frame(fmd)
coordinates(fmd0) <- c("X", "Y")
proj4string(fmd0) <- CRS("+init=epsg:27700")
fmd_sp <- as(fmd0, "SpatialPoints")
dates <- as.Date("2001-02-18")+fmd0$ReportedDay
library(spacetime)
fmd_ST <- STIDF(fmd_sp, dates, data.frame(ReportedDay=fmd0$ReportedDay))
data(SAGA_pal)
plotKML(fmd_ST, colour_scale=SAGA_pal[[1]])

## -------------- STFDF --------- ##

## results of krigeST:
library(gstat)
library(sp)
library(spacetime)
library(raster)
## define space-time variogram
sumMetricVgm <- vgmST("sumMetric",
                      space=vgm( 4.4, "Lin", 196.6,  3),
                      time =vgm( 2.2, "Lin",   1.1,  2),
                      joint=vgm(34.6, "Exp", 136.6, 12),
                      stAni=51.7)
## example from the gstat package:
data(air)
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
rr <- rural[,"2005-06-01/2005-06-03"]
rr <- as(rr,"STSDF")
x1 <- seq(from=6,to=15,by=1)
x2 <- seq(from=48,to=55,by=1)
DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
                           proj4string=CRS(proj4string(rr@sp)))
gridded(DE_gridded) <- TRUE
DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)
DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
                     modelList=sumMetricVgm)
gridded(DE_kriged@sp) <- TRUE
stplot(DE_kriged)
## plot in Google Earth:
png.width = DE_kriged@sp@grid@cells.dim[1]*20
png.height = DE_kriged@sp@grid@cells.dim[2]*20
z.lim = range(DE_kriged@data, na.rm=TRUE)
plotKML(DE_kriged, png.width=png.width, 
        png.height=png.height, z.lim=z.lim)
## add observations points:
plotKML(rr, z.lim=z.lim)

## -------------- STTDF --------- ##
library(fossil)
library(spacetime)
library(adehabitat)
data(gpxbtour)
## format the time column:
gpxbtour$ctime <- as.POSIXct(gpxbtour$time, format="%Y-%m-%dT%H:%M:%SZ")
coordinates(gpxbtour) <- ~lon+lat
proj4string(gpxbtour) <- CRS("+proj=longlat +datum=WGS84")
xy <- as.list(data.frame(t(coordinates(gpxbtour))))
gpxbtour$dist.km <- sapply(xy, function(x) { 
  deg.dist(long1=x[1], lat1=x[2], long2=xy[[1]][1], lat2=xy[[1]][2]) 
} )
## convert to a STTDF class:
gpx.ltraj <- as.ltraj(coordinates(gpxbtour), gpxbtour$ctime, id = "th")
gpx.st <- as(gpx.ltraj, "STTDF")
gpx.st$speed <- gpxbtour$speed
gpx.st@sp@proj4string <- CRS("+proj=longlat +datum=WGS84")
str(gpx.st)
plotKML(gpx.st, colour="speed")

## -------------- Spatial Metadata --------- ##
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
## subset to 20 percent:
eberg <- eberg[runif(nrow(eberg))<.1,]
eberg.md <- spMetadata(eberg["SNDMHT_A"], 
  Citation_title = 'Ebergotzen data set',
  Citation_URL = 'http://geomorphometry.org/content/ebergotzen')
plotKML(eberg["CLYMHT_A"], metadata=eberg.md)

## -------------- RasterBrickTimeSeries --------- ##
library(raster)
library(sp)
data(LST)
gridded(LST) <- ~lon+lat
proj4string(LST) <- CRS("+proj=longlat +datum=WGS84")
dates <- sapply(strsplit(names(LST), "LST"), function(x){x[[2]]})
datesf <- format(as.Date(dates, "%Y_%m_%d"), "%Y-%m-%dT%H:%M:%SZ")
## begin / end dates +/- 4 days:
TimeSpan.begin = as.POSIXct(unclass(as.POSIXct(datesf))-4*24*60*60, origin="1970-01-01") 
TimeSpan.end = as.POSIXct(unclass(as.POSIXct(datesf))+4*24*60*60, origin="1970-01-01")
## pick climatic stations in the area:
pnts <- HRtemp08[which(HRtemp08$NAME=="Pazin")[1],]
pnts <- rbind(pnts, HRtemp08[which(HRtemp08$NAME=="Crni Lug - NP Risnjak")[1],])
pnts <- rbind(pnts, HRtemp08[which(HRtemp08$NAME=="Cres")[1],])
coordinates(pnts) <- ~Lon + Lat
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
## get the dates from the file names:
LST_ll <- brick(LST[1:5])
LST_ll@title = "Time series of MODIS Land Surface Temperature images"
LST.ts <- new("RasterBrickTimeSeries", variable = "LST", sampled = pnts, 
    rasters = LST_ll, TimeSpan.begin = TimeSpan.begin[1:5], 
    TimeSpan.end = TimeSpan.end[1:5])
data(SAGA_pal)
## plot MODIS images in Google Earth:
plotKML(LST.ts, colour_scale=SAGA_pal[[1]])

## -------------- Spatial Predictions --------- ##
library(sp)
library(rgdal)
library(gstat)
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
## load grids:
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
## fit a model:
library(GSIF)
omm <- fit.gstatModel(observations = meuse, formulaString = om~dist, 
   family = gaussian(log), covariates = meuse.grid)
## produce SpatialPredictions:
om.rk <- predict(omm, predictionLocations = meuse.grid)
## plot the whole geostatical mapping project in Google Earth:
plotKML(om.rk, colour_scale = SAGA_pal[[1]], 
    png.width = gridparameters(meuse.grid)[1,"cells.dim"]*5, 
    png.height = gridparameters(meuse.grid)[2,"cells.dim"]*5)
## plot each cell as polygon:
plotKML(om.rk, colour_scale = SAGA_pal[[1]], grid2poly = TRUE)

## -------------- SpatialSamplingPattern --------- ##
library(spcosa)
library(sp)
## read a polygon map:
shpFarmsum <- readOGR(dsn = system.file("maps", package = "spcosa"), 
  layer = "farmsum")
## stratify `Farmsum' into 50 strata
myStratification <- stratify(shpFarmsum, nStrata = 50)
## sample two sampling units per stratum
mySamplingPattern <- spsample(myStratification, n = 2)
## attach the correct proj4 string:
library(RCurl)
urlExists = url.exists("http://spatialreference.org/ref/sr-org/6781/proj4/")
if(urlExists){
  nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")
  proj4string(mySamplingPattern@sample) <- CRS(nl.rd) 
  # prepare spatial domain (polygons):
  sp.domain <- as(myStratification@cells, "SpatialPolygons")
  sp.domain <- SpatialPolygonsDataFrame(sp.domain, 
     data.frame(ID=as.factor(myStratification@stratumId)), match.ID = FALSE)
  proj4string(sp.domain) <- CRS(nl.rd) 
  # create new object:
  mySamplingPattern.ssp <- new("SpatialSamplingPattern", 
     method = class(mySamplingPattern), pattern = mySamplingPattern@sample, 
     sp.domain = sp.domain)
  # the same plot now in Google Earth:
  shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
  plotKML(mySamplingPattern.ssp, shape = shape)
}

## -------------- RasterBrickSimulations --------- ##
library(sp)
library(gstat)
data(barxyz)
## define the projection system:
prj = "+proj=tmerc +lat_0=0 +lon_0=18 +k=0.9999 +x_0=6500000 +y_0=0 
  +ellps=bessel +units=m 
  +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824"
coordinates(barxyz) <- ~x+y
proj4string(barxyz) <- CRS(prj)
data(bargrid)
coordinates(bargrid) <- ~x+y
gridded(bargrid) <- TRUE
proj4string(bargrid) <- CRS(prj)
## fit a variogram and generate simulations:
Z.ovgm <- vgm(psill=1352, model="Mat", range=650, nugget=0, kappa=1.2)
sel <- runif(length(barxyz$Z))<.2  
## Note: this operation can be time consuming
sims <- krige(Z~1, barxyz[sel,], bargrid, model=Z.ovgm, nmax=20, 
   nsim=10, debug.level=-1)
## specify the cross-section:
t1 <- Line(matrix(c(bargrid@bbox[1,1], bargrid@bbox[1,2], 5073012, 5073012), ncol=2))
transect <- SpatialLines(list(Lines(list(t1), ID="t")), CRS(prj))
## glue to a RasterBrickSimulations object:
library(raster)
bardem_sims <- new("RasterBrickSimulations", variable = "elevations", 
  sampled = transect, realizations = brick(sims))
## plot the whole project and open in Google Earth:
data(R_pal)
plotKML(bardem_sims, colour_scale = R_pal[[4]])


## -------------- SpatialVectorsSimulations --------- ##
data(barstr)
data(bargrid)
library(sp)
coordinates(bargrid) <- ~ x+y
gridded(bargrid) <- TRUE
## output topology:
cell.size = bargrid@grid@cellsize[1]
bbox = bargrid@bbox
nrows = round(abs(diff(bbox[1,])/cell.size), 0) 
ncols = round(abs(diff(bbox[2,])/cell.size), 0)
gridT = GridTopology(cellcentre.offset=bbox[,1], 
  cellsize=c(cell.size,cell.size), 
  cells.dim=c(nrows, ncols))
bar_sum <- count.GridTopology(gridT, vectL=barstr[1:5])
## NOTE: this operation can be time consuming!
## plot the whole project and open in Google Earth:
plotKML(bar_sum,
    png.width = gridparameters(bargrid)[1,"cells.dim"]*5, 
    png.height = gridparameters(bargrid)[2,"cells.dim"]*5)

## -------------- SpatialMaxEntOutput --------- ##
library(maptools)
library(rgdal)
data(bigfoot)
aea.prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 
   +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
data(USAWgrids)
gridded(USAWgrids) <- ~s1+s2
proj4string(USAWgrids) <- CRS(aea.prj)
bbox <- spTransform(USAWgrids, CRS("+proj=longlat +datum=WGS84"))@bbox
sel = bigfoot$Lon > bbox[1,1] & bigfoot$Lon < bbox[1,2] &
    bigfoot$Lat > bbox[2,1] & bigfoot$Lat < bbox[2,2]
bigfoot <- bigfoot[sel,]
coordinates(bigfoot) <- ~Lon+Lat
proj4string(bigfoot) <- CRS("+proj=longlat +datum=WGS84")
library(spatstat)
bigfoot.aea <- as.ppp(spTransform(bigfoot, CRS(aea.prj)))
## Load the covariates:
sel.grids <- c("globedem","nlights03","sdroads","gcarb","twi","globcov")
library(GSIF)
library(dismo)
## run MaxEnt analysis:
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if(file.exists(jar)){
  bigfoot.smo <- MaxEnt(bigfoot.aea, USAWgrids[sel.grids])
  icon = "http://plotkml.r-forge.r-project.org/bigfoot.png"
  data(R_pal)
  plotKML(bigfoot.smo, colour_scale = R_pal[["bpy_colors"]], 
    png.width = gridparameters(USAWgrids)[1,"cells.dim"]*5, 
    png.height = gridparameters(USAWgrids)[2,"cells.dim"]*5,
    shape = icon)
}

## end of tutorial;
