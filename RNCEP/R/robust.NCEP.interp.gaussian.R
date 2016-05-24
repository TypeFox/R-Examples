robust.NCEP.interp.gaussian <-
function(variable, lat, lon, dt, 
					reanalysis2=FALSE, interpolate.space=TRUE, interpolate.time=TRUE, keep.unpacking.info=FALSE, return.units=TRUE, interp='linear', p=1){

## Latitude and longitude should be given in decimal degrees ##
## Variables must be given using the following naming conventions...
####### These variables are forecasts valid 6 hours after the reference time #########
	## 'air.2m'			==	Air Temperature				(At 2 meters)
	## 'icec.sfc'		== 	Ice Concentration				(At Surface)
	## 'pevpr.sfc'		==	Potential Evaporation Rate		(At Surface)
	## 'pres.sfc'		==	Pressure					(At Surface)
	## 'runof.sfc'		== 	Water Runoff				(At Surface)
	## 'sfcr.sfc'		==	Surface Roughness				(At Surface)
	## 'shum.2m'		==	Specific Humidity				(At 2 meters)
	## 'soilw.0-10cm'		==	Soil Moisture 				(From 0-10 cm)
	## 'soilw.10-200cm'	==	Soil Moisture 				(From 10-200 cm)
	## 'skt.sfc'		== 	Skin Temperature				(At Surface)
	## 'tmp.0-10cm'		==	Temperature of 0-10 cm layer		(From 0-10 cm)
	## 'tmp.10-200cm'		==	Temperature of 10-200 cm layer	(From 10-200 cm)
	## 'tmp.300cm'		==	Temperature at 300 cm			(From 300 cm)
	## 'uwnd.10m'		==	U-wind					(At 10 meters)
	## 'vwnd.10m'		==	V-wind					(At 10 meters)
	## 'weasd.sfc'		==	Water equivalent of snow depth	(At Surface)
########## These variables are 6 hour hindcasts from the reference time ###############
	## 'tmax.2m'		==	Maximum temperature			(At 2 meters)
	## 'tmin.2m'		==	Minimum temperature			(At 2 meters)
######### These variables are 6 hour averages starting at the reference time ##########
	## 'cfnlf.sfc'		==	Cloud forcing net longwave flux	(At Surface)
	## 'cfnsf.sfc'		==	Cloud forcing net solar flux		(At Surface)
	## 'cprat.sfc'		==	Convective precipitation rate		(At Surface)
	## 'csdlf.sfc'		==	Clear sky downward longwave flux	(At Surface)
	## 'csdsf.sfc'		==	Clear sky downward solar flux		(At Surface)
	## 'dlwrf.sfc'		==	Downward longwave radiation flux	(At Surface)
	## 'dswrf.sfc'		==	Downward solar radiation flux		(At Surface)
	## 'dswrf.ntat'		==	Downward solar radiation flux		(Nominal Top of Atmosphere)
	## 'gflux.sfc'		==	Ground heat flux				(At Surface)
	## 'lhtfl.sfc'		==	Latent heat net flux			(At Surface)
	## 'nbdsf.sfc'		==	Near IR beam downward solar flux	(At Surface)
	## 'nddsf.sfc'		==	Near IR diffuse downward solar flux	(At Surface)
	## 'nlwrs.sfc'		==	Net longwave radiation			(At Surface)
	## 'nswrs.sfc'		==	Net shortwave radiation			(At Surface)
	## 'prate.sfc'		==	Precipitation rate			(At Surface)
	## 'shtfl.sfc'		==	Sensible heat net flux			(At Surface)
	## 'uflx.sfc'		==	Momentum flux (zonal)			(At Surface)
	## 'ugwd.sfc'		==	Zonal gravity wave stress		(At Surface)
	## 'ulwrf.sfc'		==	Upward longwave radiation flux	(At Surface)
	## 'ulwrf.ntat'		==	Upward longwave radiation flux	(Nominal Top of Atmosphere)
	## 'uswrf.sfc'		==	Upward solar radiation flux		(At Surface)
	## 'uswrf.ntat'		==	Upward solar radiation flux		(Nominal Top of Atmosphere)
	## 'vbdsf.sfc'		== 	Visable beam downward solar flux	(At Surface)
	## 'vddsf.sfc'		==	Visable diffuse downward solar flux	(At Surface)
	## 'vflx.sfc'		==	Momentum flux (meridional)		(At Surface)
	## 'vgwd.sfc'		==	Meridional gravity wave stress	(At Surface)
######### These variables are 6 hour averages starting at the reference time ##########
	## 'csulf.ntat'		==	Clear Sky Upward Longwave Flux	(Nominal Top of Atmosphere)
	## 'csusf.ntat'		==	Clear Sky Upward Solar Flux		(Nominal Top of Atmosphere)
	## 'dswrf.ntat'		==	Downward Solar Radiation Flux		(Nominal Top of Atmosphere)
	## 'pres.hcb'		==	Pressure					(High Cloud Bottom)
	## 'pres.hct'		==	Pressure					(High Cloud Top)
	## 'pres.lcb'		==	Pressure					(Low Cloud Bottom)
	## 'pres.lct'		==	Pressure					(Low Cloud Top)
	## 'pres.mcb'		==	Pressure					(Middle Cloud Bottom)
	## 'pres.mct'		==	Pressure					(Middle Cloud Top)
	## 'tcdc.eatm'		==	Total Cloud Cover				(Entire Atmosphere)
	## 'ulwrf.ntat'		==	Upward Longwave Radiation Flux	(Nominal Top of Atmosphere)
	## 'uswrf.ntat'		==	Upward Solar Radiation Flux		(Nominal Top of Atmosphere)


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

#############################################################
## A lookup table for latitudes, longitudes, and variables ##
possible.lats <- c(88.542, 86.6531, 84.7532, 82.8508, 80.9473, 79.0435, 77.1394, 75.2351, 73.3307, 71.4262, 69.5217, 67.6171, 65.7125, 63.8079, 61.9033,
				59.9986, 58.0939, 56.1893, 54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, 42.8564, 40.9517, 39.047, 37.1422, 35.2375, 33.3328,
				31.4281, 29.5234, 27.6186, 25.7139, 23.8092, 21.9044, 19.9997, 18.095, 16.1902, 14.2855, 12.3808, 10.47604, 8.57131, 6.66657, 4.76184,
				2.8571, 0.952368, -0.952368, -2.8571, -4.76184, -6.66657, -8.57131, -10.47604, -12.3808, -14.2855, -16.1902, -18.095, -19.9997, -21.9044,
				-23.8092, -25.7139, -27.6186, -29.5234, -31.4281, -33.3328, -35.2375, -37.1422, -39.047, -40.9517, -42.8564, -44.7611, -46.6658, -48.5705,
				-50.4752, -52.3799, -54.2846, -56.1893, -58.0939, -59.9986, -61.9033, -63.8079, -65.7125, -67.6171, -69.5217, -71.4262, -73.3307, -75.2351,
				-77.1394, -79.0435, -80.9473, -82.8508, -84.7532, -86.6531, -88.542)
possible.lons <- seq(0, 358.125, by=1.875)
possible.variables <- c('air.2m','icec.sfc','pevpr.sfc','pres.sfc','runof.sfc','sfcr.sfc','shum.2m','soilw.0-10cm','soilw.10-200cm','skt.sfc',
				'tmp.0-10cm','tmp.10-200cm','tmp.300cm','uwnd.10m','vwnd.10m','weasd.sfc','tmax.2m','tmin.2m','cfnlf.sfc','cfnsf.sfc',
				'cprat.sfc','csdlf.sfc','csdsf.sfc','dlwrf.sfc','dswrf.sfc','dswrf.ntat','gflux.sfc','lhtfl.sfc','nbdsf.sfc','nddsf.sfc','nlwrs.sfc',
				'nswrs.sfc','prate.sfc','shtfl.sfc','uflx.sfc','ugwd.sfc','ulwrf.sfc','ulwrf.ntat','uswrf.sfc','uswrf.ntat','vbdsf.sfc','vddsf.sfc','vflx.sfc',
				'vgwd.sfc','csulf.ntat','csusf.ntat','dswrf.ntat','pres.hcb','pres.hct','pres.lcb','pres.lct','pres.mcb','pres.mct',
				'tcdc.eatm','ulwrf.ntat','uswrf.ntat')
hindcast.variables <- c('tmax.2m','tmin.2m')


######################################################
## Specify which variables are not in Reanalysis II ##
not.in.reanalysis2 <- c('sfcr.sfc','tmp.300cm','cfnlf.sfc','cfnsf.sfc','csdlf.sfc','csdlf.sfc','csdsf.sfc','nbdsf.sfc','nddsf.sfc',
	 			'nlwrs.sfc','nswrs.sfc','vbdsf.sfc', 'vddsf.sfc','csulf.ntat','csusf.ntat')

##############################################
## Specify the grid sizes in space and time ##
lat.gridsize <- abs(diff(possible.lats))*100
lon.gridsize <- 187.5
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
##### Confirm that the variable selected exists at the pressure level specified #####
if(any(possible.variables == variable[i]) == FALSE) { stop(paste("'",variable[i],"'", " not a valid variable name in reference to a gaussian grid.",sep='')) }

#######################################################################################
#### Confirm that the variable selected exists in the Reanalysis dataset specified ####
if(reanalysis2[i] == TRUE && any(not.in.reanalysis2 == variable[i])) {
	reanalysis2[i] <- FALSE
	warning(paste(variable[i]," is not in the Reanalysis II dataset.  Using Reanalysis I instead.",sep=''))
}


###################################################
## Determine what the variable 'level' should be ##
## Subdivide the variable name to extract 'name' and 'level'
name <- strsplit(variable[i], "\\.")[[1]][1]
level <- strsplit(variable[i], "\\.")[[1]][2]


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
usable.lats <- order(abs(lat[i] - possible.lats))[c(1,2)]
ilat0 <- possible.lats[max(usable.lats)]*100
ilat1 <- possible.lats[min(usable.lats)]*100
ilon0 <- floor(lon[i] * 100 / lon.gridsize) * lon.gridsize
ilon1 <- ilon0 + lon.gridsize

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
f0ts <- ifelse(variable[i] %in% hindcast.variables, 1, 0)
} else

if(interp[i] == 'linear'){
## SPACE ##
## Latitudes ##
usable.lats <- order(abs(lat[i] - possible.lats))[c(1,2)]
ilat0 <- possible.lats[max(usable.lats)]*100
ilat1 <- possible.lats[min(usable.lats)]*100
f0lat <- (lat[i] * 100.0 - ilat0) / lat.gridsize[max(usable.lats)]
f0lat <- ifelse(interpolate.space[i] == FALSE, round(f0lat, digits=0), f0lat)
## Longitudes ##
ilon0 <- floor(lon[i] * 100 / lon.gridsize) * lon.gridsize
ilon1 <- ilon0 + lon.gridsize
f0lon <- (lon[i] * 100.0 - ilon0) / lon.gridsize
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
f0ts <- ifelse(variable[i] %in% hindcast.variables, 1, 0)
} else stop("'interp' must be either 'linear' or 'IDW'")

########################################################
## Specify the area around the desired point in space ##
lat.range <- c(which(possible.lats == possible.lats[min(usable.lats)])-1, which(possible.lats == possible.lats[max(usable.lats)])-1)
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
columns <- length(seq((ilon0/100),(ilon1/100), by=lon.gridsize/100)) + 1
actual.columns <- columns-1
rows <- length(seq((ilat0/100),(ilat1/100), by=lat.gridsize[min(usable.lats)]/100))


#################################################################################
## Specify the name of the folder (depends on variable and reanalysis I or II) ##
if(reanalysis2[i] == TRUE) { gauss.name <- "gaussian_grid" } else if(any(c('2m','sfc','0-10cm','10-200cm','300cm','10m') == level)) { gauss.name <- "surface_gauss" } else { gauss.name <- "other_gauss" }


if(i == 1 | keep.unpacking.info == FALSE){
#####################################################################################################
### Query the data file online to determine the appropriate offset and scale factor to apply ##
#### FOR VARIABLES ON A GAUSSIAN GRID ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.das", sep=''), scale.offset.missingvals.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.das into a web browser to obtain an error message.", sep = ""))}
}
##
add.offset <- if(reanalysis2[i] == TRUE){as.numeric(strsplit(strsplit(grep('add_offset', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'add_offset ')[[1]][2]) } else { 0 }
scale.factor <- if(reanalysis2[i] == TRUE){as.numeric(strsplit(strsplit(grep('scale_factor', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'scale_factor ')[[1]][2]) } else { 1 }
missing.values <- as.numeric(strsplit(strsplit(grep('missing_value', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'missing_value ')[[1]][2])

var.loc.units <- min(grep(name, x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE))
all.loc.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE)
all.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = TRUE, fixed = TRUE)
units[i] <- strsplit(all.units[which(all.loc.units > var.loc.units)[1]], "\"")[[1]][2]
	
}


###################################################################################
## Download the variable for the area in space and time around the desired point ##
#### FOR VARIABLES AT A PARTICULAR PRESSURE LEVEL ####
	trying.out.west.year1 <- 1
	trying.out.west.year2 <- 1
	trying.out.east.year1 <- 1
	trying.out.east.year2 <- 1
fail <- 0
while(trying.out.west.year1 != 0){
trying.out.west.year1 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.ascii?",name,"[",beg.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[1],"]", sep=''), out.temp.west.year1), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.ascii?",name,"[",beg.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[1],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.west.year2 != 0){
trying.out.west.year2 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year2,".nc.ascii?",name,"[",end.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[1],"]", sep=''), out.temp.west.year2), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year2,".nc.ascii?",name,"[",end.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[1],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.east.year1 != 0){
trying.out.east.year1 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.ascii?",name,"[",beg.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[2],"]", sep=''), out.temp.east.year1), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year1,".nc.ascii?",name,"[",beg.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[2],"] into a web browser to obtain an error message.", sep = ""))}
}
fail <- 0
while(trying.out.east.year2 != 0){
trying.out.east.year2 <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year2,".nc.ascii?",name,"[",end.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[2],"]", sep=''), out.temp.east.year2), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2[i] == TRUE, "2",""),"/",gauss.name,"/",variable[i],".gauss.",year2,".nc.ascii?",name,"[",end.jdate,"]",ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",lat.range[2],"][",lon.range[2],"] into a web browser to obtain an error message.", sep= ""))}
}

###################################################
## Retrieve weather data from the temporary file ##
#### FOR VARIABLES ON THE GAUSSIAN GRID ####
outdata.west.year1 <- read.table(file=out.temp.west.year1, sep=',', skip=ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level),13,12), header=FALSE, na.strings=missing.values, nrows=2)
outdata.west.year2 <- read.table(file=out.temp.west.year2, sep=',', skip=ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level),13,12), header=FALSE, na.strings=missing.values, nrows=2)
outdata.east.year1 <- read.table(file=out.temp.east.year1, sep=',', skip=ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level),13,12), header=FALSE, na.strings=missing.values, nrows=2)
outdata.east.year2 <- read.table(file=out.temp.east.year2, sep=',', skip=ifelse(reanalysis2[i] == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level),13,12), header=FALSE, na.strings=missing.values, nrows=2)

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

## Calculate the standard deviation around the values ##
if(interpolate.space[i] == TRUE){
	spread[i] <- ifelse(variable[i] %in% hindcast.variables, sd(c(rec1,rec3,rec5,rec7)), sd(c(rec0,rec2,rec4,rec6)))
	} else {
	spread[i] <- NA
	}

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
outdata.west.year2 <- c()
outdata.east.year1 <- c()
outdata.east.year2 <- c()

}  ## END FOR LOOP ##

#########################################
## Disconnect from the temporary files ##
unlink(c(out.temp.west.year1, out.temp.west.year2, out.temp.east.year1, out.temp.east.year2, scale.offset.missingvals.temp))

#############################
## Return the desired data ##
return(data.frame(wx.out, units, spread))

}  ## END FUNCTION ##

