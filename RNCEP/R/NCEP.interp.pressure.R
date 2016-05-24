NCEP.interp.pressure <-
function(variable, lat, lon, dt, pressure, reanalysis2=FALSE, 
	interpolate.space=TRUE, interpolate.time=TRUE, keep.unpacking.info=FALSE,
	return.units=TRUE, interp='linear', p=1, status.bar=TRUE){

## Latitude and longitude should be given in decimal degrees ##
## 'pressure' can be one of 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10
## Variables must be given using the following naming conventions...
	## 'air'			==	Air Temperature(ºK)
	## 'hgt'			== 	Geopotential Height(m)
	## 'rhum'			==	Relative Humidity(%)
	## 'shum'			==	Specific Humidity(kg/kg)
	## 'omega' 			==	Omega [Vertical Velocity]
	## 'uwnd'			==	U-Wind Component(m/s)
	## 'vwnd'			==	V-Wind Component(m/s)


################################################################
## Determine the number of points the function will calculate ##
iterations <- max(c(length(variable),length(lat), length(lon), length(dt), length(pressure), length(reanalysis2), length(interpolate.space), length(interpolate.time), length(interp), length(p)))

## If a status bar is desired, describe the status bar parameters ##
if(status.bar){
	#importFrom(tcltk,tkProgressBar)
	#require(tcltk)
	pb <- tkProgressBar(title="Total progress", min = 0, max=iterations, width=300) 
		} else { pb <- NULL }

########################################################################################
## Recycle any variable that is shorter than the length of the longest input variable ##
variable <- rep(variable, length.out=iterations)
lat <- rep(lat, length.out=iterations)
lon <- rep(lon, length.out=iterations)
dt <- rep(dt, length.out=iterations)
pressure <- rep(pressure, length.out=iterations)
reanalysis2 <- rep(reanalysis2, length.out=iterations)
interpolate.space <- rep(interpolate.space, length.out=iterations)
interpolate.time <- rep(interpolate.time, length.out=iterations)
interp <- rep(interp, length.out=iterations)
p <- rep(p, length.out=iterations)

## If spatial interpolation is turned off (i.e. nearest neighbor interpolation), make sure that the method of interpolation is 'IDW' ##
interp <- ifelse(!interpolate.space, 'IDW', interp)

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
possible.levels <- c(1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0, 100.0, 70.0, 50.0, 30.0, 20.0, 10.0)


######################################################
## Specify which variables are not in Reanalysis II ##
not.in.reanalysis2 <- c('shum')


##############################################
## Specify the grid sizes in space and time ##
gridsize <- as.integer(250)
tgridsize <- as.integer(60*60*6) ## six hours


#######################################################
## Create temporary files to store query information ##
scale.offset.missingvals.temp <- tempfile()
out.temp <- tempfile()


#####################################################
## Create the output variable to store output data ##
wx.out <- c()
units <- c()
spread <- c()


##################################################################
## Specify that unpacking information has not yet been acquired ##
unpacking.info.acquired <- FALSE


##################################################################
## Create a function to linearly interpolate between two points ##
MK.interp <- function(pt1, pt2, pos){
	return(pt1 + pos * (pt2 - pt1))
	}


####################
####################
## Start the loop ##
for(i in 1:iterations){
##### Confirm that the variable selected exists at the pressure level specified #####
if(any(c('air','hgt','rhum','shum','omega','uwnd','vwnd') == variable[i]) == FALSE) { stop(paste("'",variable[i],"'", " not a valid variable name with reference to a pressure level",sep='')) }


##### Confirm that the specified pressure level exists #####
if(any(possible.levels == pressure[i]) == FALSE) { stop(paste(pressure[i], " not a valid pressure level.  Must be one of ", paste(possible.levels, collapse=', '),sep='')) }



#######################################################################################
#### Confirm that the variable selected exists in the Reanalysis dataset specified ####
if(reanalysis2[i] == TRUE && any(not.in.reanalysis2 == variable[i])) {
	reanalysis2[i] <- FALSE
	warning(paste(variable[i]," is not available at a specific pressure level in the Reanalysis II dataset.  Using Reanalysis I instead.",sep=''))
}


#############################################
## Convert negative longitudes to positive ##
lon[i] <- ifelse(lon[i] < 0, 360+lon[i], lon[i])


###############################################
## When interpolation requires points from both sides of the Prime Meridian, ##
## use the robust version of the function. ##
if(lon[i] >= 357.5){
	temp.out <- robust.NCEP.interp.pressure(variable=variable[i], lat=lat[i], lon=lon[i], dt=dt[i], pressure=pressure[i],
		reanalysis2=reanalysis2[i], interpolate.space=interpolate.space[i], interpolate.time=interpolate.time[i],
		keep.unpacking.info=keep.unpacking.info, return.units=return.units, interp=interp[i], p=p[i])
	wx.out[i] <- temp.out$wx.out
	spread[i] <- temp.out$spread
	if(return.units == TRUE){
		units[i] <- as.character(temp.out$units) }
	## Update the status bar ##
	if(!is.null(pb)){
		cval <- pb$getVal()
		Sys.sleep(0.000001)
		setTkProgressBar(pb, cval+1, label=paste(round((cval+1)/iterations*100, 0), "% done"))
		if(pb$getVal() == iterations) {close(pb)}
		}
	## Proceed to the next iteration ##	
	next }

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
year <- as.numeric(format(dt.f, "%Y"))
idatetime <-  floor(as.numeric(as.POSIXct(dt[i], tz='UTC')) / tgridsize)
ts0 <- idatetime * tgridsize
ts1 <- ts0 + tgridsize 
f0ts <- (as.numeric(as.POSIXct(dt[i], tz='UTC')) - ts0) / tgridsize
f0ts <- ifelse(interpolate.time[i] == FALSE, round(f0ts, digits=0), f0ts)
} else

if(interp[i] == 'linear'){
##############################################################################################
## Create some variables to detemine which points should be obtained and how to interpolate ##
## SPACE ##
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
year <- as.numeric(format(dt.f, "%Y"))
idatetime <-  floor(as.numeric(as.POSIXct(dt[i], tz='UTC')) / tgridsize)
ts0 <- idatetime * tgridsize
ts1 <- ts0 + tgridsize 
f0ts <- (as.numeric(as.POSIXct(dt[i], tz='UTC')) - ts0) / tgridsize
f0ts <- ifelse(interpolate.time[i] == FALSE, round(f0ts, digits=0), f0ts)
} else stop("'interp' must be either 'IDW' or 'linear'")

##############################################
## When interpolation requires data from two separate years, ##
## use the robust version of the function. ##
if(format(dt.f, "%m-%d %H:%M:%S") > "12-31 17:59:59") {
	temp.out <- robust.NCEP.interp.pressure(variable=variable[i], lat=lat[i], lon=lon[i], dt=dt[i], pressure=pressure[i],
		reanalysis2=reanalysis2[i], interpolate.space=interpolate.space[i], interpolate.time=interpolate.time[i],
		keep.unpacking.info=keep.unpacking.info, return.units=return.units, interp=interp[i], p=p[i])
	wx.out[i] <- temp.out$wx.out
	spread[i] <- temp.out$spread
	if(return.units == TRUE){
		units[i] <- as.character(temp.out$units) }
	next }


########################################################
## Specify the area around the desired point in space ##
lat.range <- c(which(possible.lats == (ilat1/100))-1, which(possible.lats == (ilat0/100))-1)
lon.range <- c(which(possible.lons == (ilon0/100))-1, which(possible.lons == (ilon1/100))-1)


#######################################################
## Specify the area around the desired point in time ##
## Calculate the location of the beginning and end dates in the OpenDAP file ##
beg.jdate <- as.numeric(difftime(as.POSIXct(ts0, origin='1970-01-01', tz="UTC"),
	 as.POSIXct(paste(year,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6
end.jdate <- as.numeric(difftime(as.POSIXct(ts1, origin='1970-01-01', tz="UTC"),
	 as.POSIXct(paste(year,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6


#################################################################
## Specify the expected number of columns of data for the data ##
columns <- length(seq((ilon0/100),(ilon1/100), by=2.5)) + 1
actual.columns <- columns-1
rows <- length(seq((ilat0/100),(ilat1/100), by=2.5))


if(i == 1 | keep.unpacking.info == FALSE | unpacking.info.acquired == FALSE){
#####################################################################################################
### Query the data file online to determine the appropriate offset and scale factor to apply ##
#### FOR VARIABLES AT A PARTICULAR PRESSURE LEVEL ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/pressure/",variable[i],".",year,".nc.das", sep=''), scale.offset.missingvals.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/pressure/",variable[i],".",year,".nc.das into a web browser to obtain an error message.", sep = ""))}
}
##
add.offset <- if(reanalysis2[i] == TRUE){ as.numeric(strsplit(strsplit(grep('add_offset', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'add_offset ')[[1]][2]) } else { 0 }
scale.factor <- if(reanalysis2[i] == TRUE){ as.numeric(strsplit(strsplit(grep('scale_factor', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'scale_factor ')[[1]][2]) } else { 1 }
missing.values <- as.numeric(strsplit(strsplit(grep('missing_value', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'missing_value ')[[1]][2])
if(return.units == TRUE){
	var.loc.units <- min(grep(variable[i], x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE))
	all.loc.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE)
	all.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = TRUE, fixed = TRUE)
	units[i] <- strsplit(all.units[which(all.loc.units > var.loc.units)[1]], "\"")[[1]][2]

	
}
unpacking.info.acquired <- TRUE
}


###################################################################################
## Download the variable for the area in space and time around the desired point ##
#### FOR VARIABLES AT A PARTICULAR PRESSURE LEVEL ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/pressure/",variable[i],".",year,".nc.ascii?",variable[i],"[",beg.jdate,":",end.jdate,"][",which(possible.levels == pressure[i])-1,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],":",lon.range[2],"]", sep=''), out.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/pressure/",variable[i],".",year,".nc.ascii?",variable[i],"[",beg.jdate,":",end.jdate,"][",which(possible.levels == pressure[i])-1,"][",lat.range[1],":",lat.range[2],"][",lon.range[1],":",lon.range[2],"] into a web browser to obtain an error message.", sep = ""))}
}

###################################################
## Retrieve weather data from the temporary file ##
#### FOR VARIABLES AT A PARTICULAR PRESSURE LEVEL ####
outdata <- read.table(file=out.temp, sep=',', skip=13, header=FALSE, na.strings=missing.values, nrows=((end.jdate-beg.jdate)+1)*rows)

########################################################################
## Put the weather values in the proper order, sorted by lat/lon/time ##
rec0.1 <- ifelse(outdata$V2[2] == missing.values, NA, outdata$V2[2] * scale.factor + add.offset)
rec1.1 <- ifelse(outdata$V2[4] == missing.values, NA, outdata$V2[4] * scale.factor + add.offset)
rec2.1 <- ifelse(outdata$V3[2] == missing.values, NA, outdata$V3[2] * scale.factor + add.offset)
rec3.1 <- ifelse(outdata$V3[4] == missing.values, NA, outdata$V3[4] * scale.factor + add.offset)
rec4.1 <- ifelse(outdata$V2[1] == missing.values, NA, outdata$V2[1] * scale.factor + add.offset)
rec5.1 <- ifelse(outdata$V2[3] == missing.values, NA, outdata$V2[3] * scale.factor + add.offset)
rec6.1 <- ifelse(outdata$V3[1] == missing.values, NA, outdata$V3[1] * scale.factor + add.offset)
rec7.1 <- ifelse(outdata$V3[3] == missing.values, NA, outdata$V3[3] * scale.factor + add.offset)

######################################
## Interpolate the weather variable ##
if(interp[i] == 'IDW' | interp[i] == 'idw'){
t1 <- sum((lat0.lon0.f*rec0.1) + (lat0.lon1.f*rec2.1) + (lat1.lon0.f*rec4.1) + (lat1.lon1.f*rec6.1))
t2 <- sum((lat0.lon0.f*rec1.1) + (lat0.lon1.f*rec3.1) + (lat1.lon0.f*rec5.1) + (lat1.lon1.f*rec7.1))
wx.out[i] <- MK.interp(t1, t2, f0ts)

} else {

## First in latitude ##
rec0 <- MK.interp(rec0.1, rec4.1, f0lat)
rec1 <- MK.interp(rec1.1, rec5.1, f0lat)
rec2 <- MK.interp(rec2.1, rec6.1, f0lat)
rec3 <- MK.interp(rec3.1, rec7.1, f0lat)
## Then in longitude ##
rec0 <- MK.interp(rec0, rec2, f0lon)
rec1 <- MK.interp(rec1, rec3, f0lon)
## Finally in time ##
wx.out[i] <- MK.interp(rec0, rec1, f0ts)

}

## Calculate the standard deviation of the values ##
if(interpolate.space[i] == TRUE){
	if(interpolate.time[i] == TRUE){
		spread[i] <- sd(c(rec0.1,rec1.1,rec2.1,rec3.1,rec4.1,rec5.1,rec6.1,rec7.1))
		} else 
	if(f0ts == 1){
		spread[i] <- sd(c(rec1.1,rec3.1,rec5.1,rec7.1))
		} else 
		spread[i] <- sd(c(rec0.1,rec2.1,rec4.1,rec6.1))
	} else 
if(interpolate.space[i] == FALSE){
	if(interpolate.time[i] == TRUE){
		spread[i] <- sd(c(t1, t2))
		} else 
	if(interpolate.time[i] == FALSE){
		spread[i] <- NA
		}
	}

#################################################
## Clear some variables for the next iteration ##
rec0.1 <- c()
rec1.1 <- c()
rec2.1 <- c()
rec3.1 <- c()
rec4.1 <- c()
rec5.1 <- c()
rec6.1 <- c()
rec7.1 <- c()
outdata <- c()

## Update the status bar ##
	if(!is.null(pb)){
		cval <- pb$getVal()
		Sys.sleep(0.000001)
		setTkProgressBar(pb, cval+1, label=paste(round((cval+1)/iterations*100, 0), "% done"))
		}


}  ## END FOR LOOP ##

#########################################
## Disconnect from the temporary files ##
unlink(c(scale.offset.missingvals.temp, out.temp))

##########################
## Close the status bar ##
if(!is.null(pb)) { if(pb$getVal() == iterations) {close(pb)} }

#####################
## Print the units ##
if(return.units == TRUE){
units <- units[is.na(units)==FALSE]
	for(x in 1:length(unique(units))){
		print(noquote(paste("Units of variable '", unique(variable)[x], "' are ", unique(units)[x], sep='')))
	}
}

#############################
## Return the desired data ##
attr(wx.out, "standard deviation") <- spread
return(wx.out)

}  ## END FUNCTION ##

