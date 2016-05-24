robust.NCEP.interp.surface <-
function(variable, lat, lon, dt, 
					reanalysis2=FALSE, interpolate.space=TRUE, interpolate.time=TRUE, keep.unpacking.info=FALSE, return.units=TRUE, interp='linear', p=1){

## Latitude and longitude should be given in decimal degrees ##
## Variables must be given using the following naming conventions...
	## 'air.sig995'		==	Air Temperature				(Near Surface)
	## 'lftx.sfc'		==	Surface Lifted Index (ºK)		(At Surface)
	## 'lftx4.sfc' 		==	Best (4-layer) Lifted Index(ºK)	(At Surface)
	## 'omega.sig995'		==	Omega [Vertical Velocity]		(Near Surface)
	## 'pottmp.sig995'	==	Potential Temperature(ºK)		(Near Surface)
	## 'pr_wtr.eatm'		==	Precipitable Water(kg/m2)		(Entire Atmosphere)
	## 'pres.sfc'		==	Pressure					(At Surface)
	## 'rhum.sig995'		==	Relative Humidity				(Near Surface)
	## 'slp'			==	Sea Level Pressure(Pa)			(Sea Level)
	## 'mslp'			==	Mean Sea Level Pressure			(Sea Level)
	## 'uwnd.sig995'		==	U-Wind Component (m/s)			(Near Surface)
	## 'vwnd.sig995'		==	V-Wind Component (m/s)			(Near Surface)


################################################################
## Determine the number of points the function will calculate ##
iterations <- max(c(length(variable),length(lat), length(lon), length(dt), length(reanalysis2), length(interpolate.space), length(interpolate.time), length(interp), length(p)))


########################################################################################
## Recycle any variable that is shorter than the length of the longest input variable ##
variable <- rep(variable, length.out=iterations)
lat <- rep(lat, length.out=iterations)
lon <- rep(lon, length.out=iterations)
dt <- rep(dt, length.out=iterations)
reanalysis2 <- rep(reanalysis2, length.out=iterations)
interpolate.space <- rep(interpolate.space, length.out=iterations)
interpolate.time <- rep(interpolate.time, length.out=iterations)
interp <- rep(interp, length.out=iterations)
p <- rep(p, length.out=iterations)

#################################################
## A lookup table for latitudes and longitudes ##
###################################################################
## A lookup table for latitudes, longitudes, and pressure levels ##
possible.lats <- c(90.0, 87.5, 85.0, 82.5, 80.0, 77.5, 75.0, 72.5, 70.0, 67.5, 65.0, 62.5, 60.0, 57.5, 55.0, 52.5, 50.0, 47.5, 45.0, 42.5, 40.0, 37.5, 35.0, 32.5, 30.0, 27.5, 25.0, 22.5, 20.0, 17.5, 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0.0, -2.5, -5.0, -7.5, -10.0, -12.5, -15.0, -17.5, -20.0, -22.5, -25.0, -27.5, -30.0, -32.5, -35.0, -37.5, -40.0, -42.5, -45.0, -47.5, -50.0, -52.5, -55.0, -57.5, -60.0, -62.5, -65.0, -67.5, -70.0, -72.5, -75.0, -77.5, -80.0, -82.5, -85.0, -87.5, -90.0)
possible.lons <- c(0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0,
				47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 
				92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5, 110.0, 112.5, 115.0, 117.5, 120.0, 122.5, 125.0, 127.5, 130.0, 
				132.5, 135.0, 137.5, 140.0, 142.5, 145.0, 147.5, 150.0, 152.5, 155.0, 157.5, 160.0, 162.5, 165.0, 167.5, 170.0, 
				172.5, 175.0, 177.5, 180.0, 182.5, 185.0, 187.5, 190.0, 192.5, 195.0, 197.5, 200.0, 202.5, 205.0, 207.5, 210.0, 
				212.5, 215.0, 217.5, 220.0, 222.5, 225.0, 227.5, 230.0, 232.5, 235.0, 237.5, 240.0, 242.5, 245.0, 247.5, 250.0, 
				252.5, 255.0, 257.5, 260.0, 262.5, 265.0, 267.5, 270.0, 272.5, 275.0, 277.5, 280.0, 282.5, 285.0, 287.5, 290.0, 
				292.5, 295.0, 297.5, 300.0, 302.5, 305.0, 307.5, 310.0, 312.5, 315.0, 317.5, 320.0, 322.5, 325.0, 327.5, 330.0, 
				332.5, 335.0, 337.5, 340.0, 342.5, 345.0, 347.5, 350.0, 352.5, 355.0, 357.5)


######################################################
## Specify which variables are not in Reanalysis II ##
not.in.reanalysis2 <- c('air.sig995','lftx.sfc','lftx4.sfc','omega.sig995','pottmp.sig995','rhum.sig995','slp','uwnd.sig995','vwnd.sig995')
not.in.reanalysis1 <- c('mslp')

##############################################
## Specify the grid sizes in space and time ##
gridsize <- as.integer(250)
tgridsize <- as.integer(60*60*6) ## six hours


#######################################################
## Create temporary files to store query information ##
scale.offset.missingvals.temp <- tempfile()
out.temp.west.year1 <- tempfile()
out.temp.west.year2 <- tempfile()
out.temp.east.year1 <- tempfile()
out.temp.east.year2 <- tempfile()


#####################################################
## Create the output variable to store output data ##
wx.out <- c()
units <- c()
spread <- c()

##################################################################
## Create a function to linearly interpolate between two points ##
MK.interp <- function(pt1, pt2, pos){
	return(pt1 + pos * (pt2 - pt1))
	}


####################
####################
## Start the loop ##
for(i in 1:iterations){
##### Confirm that the variable selected exists at the level specified #####
if(any(c('air.sig995','lftx.sfc','lftx4.sfc','omega.sig995','pottmp.sig995','pr_wtr.eatm','pres.sfc','rhum.sig995','slp','mslp','uwnd.sig995','vwnd.sig995') == variable[i]) == FALSE) { stop(paste("'",variable[i],"'", " not a valid variable name in reference to the surface.",sep='')) }


#######################################################################################
#### Confirm that the variable selected exists in the Reanalysis dataset specified ####
if(reanalysis2[i] == TRUE && any(not.in.reanalysis2 == variable[i])) {
	reanalysis2[i] <- FALSE
	warning(paste(variable[i]," is not available at the surface in the Reanalysis II dataset.  Using Reanalysis I instead.",sep=''))
}
if(reanalysis2[i] == FALSE && any(not.in.reanalysis1 == variable[i])) {
	reanalysis2[i] <- TRUE
	warning(paste(variable[i]," is not available at the surface in the Reanalysis I dataset.  Using Reanalysis II instead.",sep=''))
}


###################################################
## Determine what the variable 'level' should be ##
## Subdivide the variable name to extract 'name' and 'level'
name <- strsplit(variable[i], "\\.")[[1]][1]
level <- ifelse(length(strsplit(variable[i], "\\.")[[1]]) > 1, strsplit(variable[i], "\\.")[[1]][2], "")
	

#############################################
## Convert negative longitudes to positive ##
lon[i] <- ifelse(lon[i] < 0, 360+lon[i], lon[i])


#######################################################################################################
##### Determine if the request will need data from both sides of the prime meridian to interpolate ####
if(lon[i] >= 357.5){
	cross.prime <- TRUE
	} else { cross.prime <- FALSE }


###########################################
if(interp[i] == 'IDW' | interp[i] == 'idw'){
NCEP.deg.dist <- function (long1, lat1, long2, lat2) 
{
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- long1 * rad
    b1 <- lat2 * rad
    b2 <- long2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 40041.47/(2 * pi)
    d <- R * c
    return(d)
}

##############################################################################################
## Create some variables to detemine which points should be obtained and how to interpolate ##
## SPACE ##
## Determine the gridpoints surrounding the desired location ##
ilat0 <- floor(lat[i] * 100 / gridsize) * gridsize
ilat1 <- ilat0 + gridsize
ilon0 <- floor(lon[i] * 100 / gridsize) * gridsize
ilon1 <- ilon0 + gridsize

## Calculate the distance of each gridpoint from the desired location on a great circle ##
lat0.lon0.d <- NCEP.deg.dist(long1=ilon0/100, lat1=ilat0/100, long2=lon[i], lat2=lat[i])^p[i]
lat0.lon1.d <- NCEP.deg.dist(long1=ilon1/100, lat1=ilat0/100, long2=lon[i], lat2=lat[i])^p[i]
lat1.lon0.d <- NCEP.deg.dist(long1=ilon0/100, lat1=ilat1/100, long2=lon[i], lat2=lat[i])^p[i]
lat1.lon1.d <- NCEP.deg.dist(long1=ilon1/100, lat1=ilat1/100, long2=lon[i], lat2=lat[i])^p[i]

## Make sure that the interpolated point doesn't fall on an existing grid point ##
if(any(c(lat0.lon0.d, lat0.lon1.d, lat1.lon0.d, lat1.lon1.d) == 0)){
	min.d <- which(c(lat0.lon0.d, lat0.lon1.d, lat1.lon0.d, lat1.lon1.d) == 0)
	interpolate.space[i] <- FALSE
	} else {
## Take the inverse of those distances ##
lat0.lon0.inv.d <- 1/lat0.lon0.d
lat0.lon1.inv.d <- 1/lat0.lon1.d
lat1.lon0.inv.d <- 1/lat1.lon0.d
lat1.lon1.inv.d <- 1/lat1.lon1.d

## Calculate the total of all of the inverse distances ##
total.inv.d <- sum(lat0.lon0.inv.d, lat0.lon1.inv.d, lat1.lon0.inv.d, lat1.lon1.inv.d)

## Determine which point is closest if nearest neighbor interpolation is desired ##
all.latlons <- c(lat0.lon0.inv.d, lat0.lon1.inv.d, lat1.lon0.inv.d, lat1.lon1.inv.d)
min.d <- which(all.latlons == max(all.latlons))
	}

## Determine the relative influence of each gridpoint based on its distance from the desired location ##
lat0.lon0.f <- if(interpolate.space[i] == TRUE) { lat0.lon0.inv.d/total.inv.d } else if(min.d == 1) {1} else {0} 
lat0.lon1.f <- if(interpolate.space[i] == TRUE) { lat0.lon1.inv.d/total.inv.d } else if(min.d == 2) {1} else {0}
lat1.lon0.f <- if(interpolate.space[i] == TRUE) { lat1.lon0.inv.d/total.inv.d } else if(min.d == 3) {1} else {0}
lat1.lon1.f <- if(interpolate.space[i] == TRUE) { lat1.lon1.inv.d/total.inv.d } else if(min.d == 4) {1} else {0}

## TIME ##
dt.f <- strptime(dt[i], "%Y-%m-%d %H:%M:%S",'UTC')
###########################################################################
## Determine whether two years of data will be needed for interpolation  ##
if(format(dt.f, "%m-%d %H:%M:%S") > "12-31 17:59:59") {
	two.years <- TRUE
	} else {two.years <- FALSE}
year1 <- as.numeric(format(dt.f, "%Y"))
year2 <- as.numeric(format(dt.f, "%Y")) + ifelse(two.years == TRUE, 1, 0)
idatetime <-  floor(as.numeric(as.POSIXct(dt[i], tz='UTC')) / tgridsize)
ts0 <- idatetime * tgridsize
ts1 <- ts0 + tgridsize 
f0ts <- (as.numeric(as.POSIXct(dt[i], tz='UTC')) - ts0) / tgridsize
f0ts <- ifelse(interpolate.time[i] == FALSE, round(f0ts, digits=0), f0ts)
} else

if(interp[i] == 'linear'){
## Latitudes ##
ilat0 <- floor(lat[i] * 100 / gridsize) * gridsize
ilat1 <- ilat0 + gridsize
f0lat <- (lat[i] * 100.0 - ilat0) / gridsize
f0lat <- ifelse(interpolate.space[i] == FALSE, round(f0lat, digits=0), f0lat)
## Longitudes ##
ilon0 <- floor(lon[i] * 100 / gridsize) * gridsize
ilon1 <- ilon0 + gridsize
f0lon <- (lon[i] * 100.0 - ilon0) / gridsize;
f0lon <- ifelse(interpolate.space[i] == FALSE, round(f0lon, digits=0), f0lon)

## TIME ##
dt.f <- strptime(dt[i], "%Y-%m-%d %H:%M:%S",'UTC')
###########################################################################
## Determine whether two years of data will be needed for interpolation  ##
if(format(dt.f, "%m-%d %H:%M:%S") > "12-31 17:59:59") {
	two.years <- TRUE
	} else {two.years <- FALSE}
year1 <- as.numeric(format(dt.f, "%Y"))
year2 <- as.numeric(format(dt.f, "%Y")) + ifelse(two.years == TRUE, 1, 0)
idatetime <-  floor(as.numeric(as.POSIXct(dt[i], tz='UTC')) / tgridsize)
ts0 <- idatetime * tgridsize
ts1 <- ts0 + tgridsize 
f0ts <- (as.numeric(as.POSIXct(dt[i], tz='UTC')) - ts0) / tgridsize
f0ts <- ifelse(interpolate.time[i] == FALSE, round(f0ts, digits=0), f0ts)
} else stop("'interp' must be either 'IDW' or 'linear'")

########################################################
## Specify the area around the desired point in space ##
lat.range <- c(which(possible.lats == (ilat1/100))-1, which(possible.lats == (ilat0/100))-1)
lon.range <- if(cross.prime == TRUE){c(143, 0)} else {c(which(possible.lons == (ilon0/100))-1, which(possible.lons == (ilon1/100))-1)}


#######################################################
## Specify the area around the desired point in time ##
## Calculate the location of the beginning and end dates in the OpenDAP file ##
beg.jdate <- as.numeric(difftime(as.POSIXct(ts0, origin='1970-01-01', tz="UTC"),
	 as.POSIXct(paste(year1,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6
end.jdate <- as.numeric(difftime(as.POSIXct(ts1, origin='1970-01-01', tz="UTC"),
	 as.POSIXct(paste(year2,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6


#################################################################
## Specify the expected number of columns of data for the data ##
columns <- length(seq((ilon0/100),(ilon1/100), by=2.5)) + 1
actual.columns <- columns-1
rows <- length(seq((ilat0/100),(ilat1/100), by=2.5))


if(i == 1 | keep.unpacking.info == FALSE){
#####################################################################################################
### Query the data file online to determine the appropriate offset and scale factor to apply ##
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.das", sep=''), scale.offset.missingvals.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.das into a web browser to obtain an error message.", sep = ""))}
}
##
add.offset <- if(reanalysis2[i] == TRUE){ as.numeric(strsplit(strsplit(grep('add_offset', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'add_offset ')[[1]][2]) } else { 0 }
scale.factor <- if(reanalysis2[i] == TRUE){ as.numeric(strsplit(strsplit(grep('scale_factor', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'scale_factor ')[[1]][2]) } else { 1 }
missing.values <- as.numeric(strsplit(strsplit(grep('missing_value', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'missing_value ')[[1]][2])

var.loc.units <- min(grep(name, x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE))
all.loc.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE)
all.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = TRUE, fixed = TRUE)
units[i] <- strsplit(all.units[which(all.loc.units > var.loc.units)[1]], "\"")[[1]][2]

}


###################################################################################
## Download the variable for the area in space and time around the desired point ##
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
	trying.out.west.year1 <- 1
	trying.out.west.year2 <- 1
	trying.out.east.year1 <- 1
	trying.out.east.year2 <- 1
fail <- 0
while(trying.out.west.year1 != 0){
trying.out.west.year1 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.ascii?",name,"[",beg.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],"]", sep=''), out.temp.west.year1), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.ascii?",name,"[",beg.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.west.year2 != 0){
trying.out.west.year2 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year2,".nc.ascii?",name,"[",end.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],"]", sep=''), out.temp.west.year2), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year2,".nc.ascii?",name,"[",end.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.east.year1 != 0){
trying.out.east.year1 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.ascii?",name,"[",beg.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[2],"]", sep=''), out.temp.east.year1), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\n Try entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year1,".nc.ascii?",name,"[",beg.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[2],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.east.year2 != 0){
trying.out.east.year2 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year2,".nc.ascii?",name,"[",end.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[2],"]", sep=''), out.temp.east.year2), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/surface/",variable[i],".",year2,".nc.ascii?",name,"[",end.jdate,"][",lat.range[1],":",lat.range[2],"][",lon.range[2],"] into a web browser to obtain an error message.", sep = ""))}
}


####################################################
## Retrieve weather data from the temporary files ##
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
outdata.west.year1 <- read.table(file=out.temp.west.year1, sep=',', skip=12, header=FALSE, na.strings=missing.values, nrows=2)
outdata.west.year2 <- read.table(file=out.temp.west.year2, sep=',', skip=12, header=FALSE, na.strings=missing.values, nrows=2)
outdata.east.year1 <- read.table(file=out.temp.east.year1, sep=',', skip=12, header=FALSE, na.strings=missing.values, nrows=2)
outdata.east.year2 <- read.table(file=out.temp.east.year2, sep=',', skip=12, header=FALSE, na.strings=missing.values, nrows=2)


########################################################################
## Put the weather values in the proper order, sorted by lat/lon/time ##
rec0 <- ifelse(outdata.west.year1$V2[2] == missing.values, NA, outdata.west.year1$V2[2] * scale.factor + add.offset)
rec1 <- ifelse(outdata.west.year2$V2[2] == missing.values, NA, outdata.west.year2$V2[2] * scale.factor + add.offset)
rec2 <- ifelse(outdata.east.year1$V2[2] == missing.values, NA, outdata.east.year1$V2[2] * scale.factor + add.offset)
rec3 <- ifelse(outdata.east.year2$V2[2] == missing.values, NA, outdata.east.year2$V2[2] * scale.factor + add.offset)
rec4 <- ifelse(outdata.west.year1$V2[1] == missing.values, NA, outdata.west.year1$V2[1] * scale.factor + add.offset)
rec5 <- ifelse(outdata.west.year2$V2[1] == missing.values, NA, outdata.west.year2$V2[1] * scale.factor + add.offset)
rec6 <- ifelse(outdata.east.year1$V2[1] == missing.values, NA, outdata.east.year1$V2[1] * scale.factor + add.offset)
rec7 <- ifelse(outdata.east.year2$V2[1] == missing.values, NA, outdata.east.year2$V2[1] * scale.factor + add.offset)

######################################
## Interpolate the weather variable ##
if(interp[i] == 'IDW' | interp[i] == 'idw'){
t1 <- sum((lat0.lon0.f*rec0) + (lat0.lon1.f*rec2) + (lat1.lon0.f*rec4) + (lat1.lon1.f*rec6))
t2 <- sum((lat0.lon0.f*rec1) + (lat0.lon1.f*rec3) + (lat1.lon0.f*rec5) + (lat1.lon1.f*rec7))
wx.out[i] <- MK.interp(t1, t2, f0ts)

} else {

## First in latitude ##
rec0 <- MK.interp(rec0, rec4, f0lat)
rec1 <- MK.interp(rec1, rec5, f0lat)
rec2 <- MK.interp(rec2, rec6, f0lat)
rec3 <- MK.interp(rec3, rec7, f0lat)
## Then in longitude ##
rec0 <- MK.interp(rec0, rec2, f0lon)
rec1 <- MK.interp(rec1, rec3, f0lon)
## Finally in time ##
wx.out[i] <- MK.interp(rec0, rec1, f0ts)

## Calculate the standard deviation of the values ##
if(interpolate.space[i] == TRUE){
	if(interpolate.time[i] == TRUE){
		spread[i] <- sd(c(outdata.west.year1$V2[c(1,2)],outdata.west.year2$V2[c(1,2)],outdata.east.year1$V2[c(1,2)],outdata.east.year2$V2[c(1,2)]))
		} else 
	if(f0ts == 1){
		spread[i] <- sd(c(outdata.west.year2$V2[c(1,2)],outdata.east.year2$V2[c(1,2)]))
		} else
		spread[i] <- sd(c(outdata.west.year1$V2[c(1,2)],outdata.east.year1$V2[c(1,2)]))
	} else 
if(interpolate.space[i] == FALSE){
	if(interpolate.time[i] == TRUE){
		spread[i] <- sd(c(t1,t2)) 
	} else {
		spread[i] <- NA
		}
}
}

#################################################
## Clear some variables for the next iteration ##
rec0 <- c()
rec1 <- c()
rec2 <- c()
rec3 <- c()
rec4 <- c()
rec5 <- c()
rec6 <- c()
rec7 <- c()
outdata.west.year1 <- c()
outdata.east.year1 <- c()
outdata.west.year2 <- c()
outdata.east.year2 <- c()


}  ## END FOR LOOP ##

#########################################
## Disconnect from the temporary files ##
unlink(c(out.temp.west.year1, out.temp.west.year2, out.temp.east.year1, out.temp.east.year2, scale.offset.missingvals.temp))


#############################
## Return the desired data ##
return(data.frame(wx.out,units,spread))

}  ## END FUNCTION ##

