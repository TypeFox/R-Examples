NCEP.gather.gaussian <-
function(variable ,months.minmax, years.minmax, lat.minmax,
	lon.minmax,  reanalysis2=FALSE, return.units=TRUE, increments=NULL, pb=NULL) {


## Latitude and longitude should be given in decimal degrees ##
## 'months.minmax' must be given in the numeric format ##
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

######################################################
## Specify which variables are not in Reanalysis II ##
not.in.reanalysis2 <- c('sfcr.sfc','tmp.300cm','cfnlf.sfc','cfnsf.sfc','csdlf.sfc','csdlf.sfc','csdsf.sfc','nbdsf.sfc','nddsf.sfc',
	 			'nlwrs.sfc','nswrs.sfc','vbdsf.sfc', 'vddsf.sfc','csulf.ntat','csusf.ntat')


##########################################################################
##### Confirm that the variable selected exists on the gaussian grid #####
if(any(possible.variables == variable) == FALSE) { stop(paste("'",variable,"'", " not a valid variable name in reference to a gaussian grid.",sep='')) }


#######################################################################################
#### Confirm that the variable selected exists in the Reanalysis dataset specified ####
if(reanalysis2 == TRUE && any(not.in.reanalysis2 == variable)) {
	reanalysis2 <- FALSE
	warning("This variable is not available on the gaussian grid in the Reanalysis II dataset.  Using Reanalysis I instead.")
}

## Figure out which latitudes and longitudes should be obtained such that the lat and lon input from the user is retained ##
## Latitudes ##
if(any(possible.lats == lat.minmax[1]) == FALSE) { 
usable.lats.min <- order(abs(lat.minmax[1] - possible.lats))[c(1,2)]
lat.minmax[1] <- possible.lats[max(usable.lats.min)]
}
if(any(possible.lats == lat.minmax[2]) == FALSE) {
usable.lats.max <- order(abs(lat.minmax[2] - possible.lats))[c(1,2)]
lat.minmax[2] <- possible.lats[min(usable.lats.max)]
}
## Longitude ##
lon.minmax[1] <- floor(lon.minmax[1] / 1.875) * 1.875
lon.minmax[2] <- ifelse((ceiling(lon.minmax[2] / 1.875) * 1.875) > 358.125, 358.125, (ceiling(lon.minmax[2] / 1.875) * 1.875))


#############################
## Now determine which points should be included from specified location ##
lat.range <- c(which(possible.lats == lat.minmax[2])-1, which(possible.lats == lat.minmax[1])-1)
lon.range <- c(which(possible.lons == lon.minmax[1])-1, which(possible.lons == lon.minmax[2])-1)


#####################################
## Determine the beginning and end days for the months and years ##
years <- seq(years.minmax[1],years.minmax[length(years.minmax)])
months <- seq(months.minmax[1], months.minmax[length(months.minmax)])
end.day <- c()
order <- 1
for(i in years){
for(m in months){
if (as.numeric(format(as.Date(paste('1',m,i, sep='/'), format('%d/%m/%Y'), tz='UTC') + 30, '%m')) == m){
	end.day[order] <- 31 
	order <- order + 1 } else
if (as.numeric(format(as.Date(paste('1',m,i, sep='/'), format('%d/%m/%Y'), tz='UTC') + 29, '%m')) == m){
	end.day[order] <- 30 
	order <- order + 1  } else
if (as.numeric(format(as.Date(paste('1',m,i, sep='/'), format('%d/%m/%Y'), tz='UTC') + 28, '%m')) == m){
	end.day[order] <- 29 
	order <- order + 1  } else
if (as.numeric(format(as.Date(paste('1',m,i, sep='/'), format('%d/%m/%Y'), tz='UTC') + 27, '%m')) == m){
	end.day[order] <- 28 
	order <- order + 1  }
}}
beg.day <- rep(1, length(end.day))


##########################################
## Create a vector of names to represent the times involved (year.month.day.hour)
year.names <- c()
loop <- 1
for(y in 1:length(years)){
year.names <- append(year.names, rep(years[y], sum(end.day[loop:(loop+length(months)-1)])*4))
loop <- loop + length(months)
}
##
month.names <- c()
loop <- 1
for(y in 1:length(years)){
for(m in 1:length(months)){
month.names <- append(month.names, sprintf("%02d", rep(months[m], each=end.day[loop]*4)))
loop <- loop + 1
}}
##
day.names <- c()
for(i in 1:length(end.day)){
	day.names <- append(day.names, sprintf("%02d", rep(seq(1,end.day[i]), each=4)))
	
}
hour.names <- sprintf("%02d", rep(seq(0,18, length.out=4), length(day.names)/4))

time.names <- paste(year.names, month.names, day.names, hour.names, sep='_')

		
###########################################
## Create an empty matrix using the dimensions of the x.coord and y.coord variables to store the output data ##
out.wx.data <- array(data=NA, dim=c(length(possible.lats[which(possible.lats == lat.minmax[2]):which(possible.lats == lat.minmax[1])]), length(possible.lons[which(possible.lons == lon.minmax[1]):which(possible.lons == lon.minmax[2])]), length(time.names)), 
	dimnames=list(possible.lats[which(possible.lats == lat.minmax[2]):which(possible.lats == lat.minmax[1])],possible.lons[which(possible.lons == lon.minmax[1]):which(possible.lons == lon.minmax[2])], time.names))


##########################################
## Create an empty vector ##
outdata <- c()


##########################################
## Create temporary files to store query information ##
scale.offset.missingvals.temp <- tempfile()
out.temp <- tempfile()


###################################################
## Determine what the variable 'level' should be ##
## Subdivide the variable name to extract 'name' and 'level'
name <- strsplit(variable, "\\.")[[1]][1]
level <- ifelse(length(strsplit(variable, "\\.")[[1]]) > 1, strsplit(variable, "\\.")[[1]][2], "")

#################################################################################
## Specify the name of the folder (depends on variable and reanalysis I or II) ##
if(reanalysis2 == TRUE) { gauss.name <- "gaussian_grid" } else if(any(c('2m','sfc','0-10cm','10-200cm','300cm','10m') == level)) { gauss.name <- "surface_gauss" } else { gauss.name <- "other_gauss" }


#################################
### Begin the Loop ##############
loop.num <- 1
observations <- c()
for(year in years){

## Calculate the location of the beginning and end dates in the OpenDAP file ##
beg.jdate <- as.numeric(difftime(as.POSIXct(paste(year,"/",months[1],"/",beg.day[loop.num]," ","0:0:0", sep=''), format="%Y/%m/%d %H:%M:%S", tz="UTC"),
	 as.POSIXct(paste(year,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6
end.jdate <- as.numeric(difftime(as.POSIXct(paste(year,"/",months[length(months)],"/",end.day[loop.num*length(months)]," ","18:0:0", sep=''), format="%Y/%m/%d %H:%M:%S", tz="UTC"),
	 as.POSIXct(paste(year,"/1/1 0:0:0", sep=''), "%Y/%m/%d %H:%M:%S", tz="UTC"), units='hours'))/6

## Specify the expected number of columns of data for each year ##
columns <- length(possible.lons[which(possible.lons == lon.minmax[1]):which(possible.lons == lon.minmax[2])]) + 1
actual.columns <- columns-1
rows <- length(possible.lats[which(possible.lats == lat.minmax[2]):which(possible.lats == lat.minmax[1])])


if(loop.num == 1){
### First query the data file online to determine the appropriate offset and scale factor to apply ##
#### FOR VARIABLES ON THE GAUSSIAN GRID ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/",gauss.name,"/",variable,".gauss.",year,".nc.das", sep=''), scale.offset.missingvals.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided. 
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/",gauss.name,"/",variable,".gauss.",year,".nc.das into a web browser to obtain an error message.",sep=""))}
}

add.offset <- if(reanalysis2 == TRUE){ as.numeric(strsplit(strsplit(grep('add_offset', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'add_offset ')[[1]][2]) } else { 0 }
scale.factor <- if(reanalysis2 == TRUE){ as.numeric(strsplit(strsplit(grep('scale_factor', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'scale_factor ')[[1]][2]) } else { 1 }
missing.values <- as.numeric(strsplit(strsplit(grep('missing_value', x=readLines(scale.offset.missingvals.temp), value=TRUE, fixed=TRUE), ';')[[1]][1], 'missing_value ')[[1]][2])
if(return.units == TRUE){
	var.loc.units <- min(grep(name, x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE))
	all.loc.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = FALSE, fixed = TRUE)
	all.units <- grep("String units", x = readLines(scale.offset.missingvals.temp), value = TRUE, fixed = TRUE)
	units <- strsplit(all.units[which(all.loc.units > var.loc.units)[1]], "\"")[[1]][2]

}
}


## Download the variable for the desired year, month(s), and location ##
#### FOR VARIABLES ON THE GAUSSIAN GRID ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2", ""),"/",gauss.name,"/",variable,".gauss.",year,".nc.ascii?",name,"[",beg.jdate,":",end.jdate,"]",ifelse(reanalysis2 == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",ifelse(length(lat.range) > 1, lat.range[2], lat.range[1]),"][",lon.range[1],":",ifelse(length(lon.range) > 1, lon.range[2], lon.range[1]),"]", sep=''), out.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2", ""),"/",gauss.name,"/",variable,".gauss.",year,".nc.ascii?",name,"[",beg.jdate,":",end.jdate,"]",ifelse(reanalysis2 == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level), "[0]",""),"[",lat.range[1],":",ifelse(length(lat.range) > 1, lat.range[2], lat.range[1]),"][",lon.range[1],":",ifelse(length(lon.range) > 1, lon.range[2], lon.range[1]),"] into a web browser to obtain an error message.",sep=""))}
}


## Retrieve weather data from the temporary file ##
## FOR VARIABLES ON THE GAUSSIAN GRID ##
outdata <- read.table(file=out.temp, sep=',', skip=ifelse(reanalysis2 == TRUE && all(c('sfc','ntat','hcb','hct','lcb','lct','mcb','mct','eatm') != level),13,12), header=FALSE, na.strings=missing.values, nrows=((end.jdate-beg.jdate)+1)*rows)


## Determine the number of observations (time)
observations[loop.num] <- ifelse(loop.num == 1, length(outdata$V2)/rows, observations[loop.num-1] + length(outdata$V2)/rows)


## Determine the number of observations for the individual year ##
## And go through each one individually ##
obs <- seq(ifelse(loop.num==1,1,observations[loop.num-1]+1),observations[loop.num])

for (i in 1:length(obs)){
n <- seq(1,length(outdata$V2), by=rows)

t.out <- outdata[c(seq(n[i],n[i]+rows-1)),c(2:columns)] * scale.factor + add.offset

## Set all of the missing values to NA ##
t.out[t.out == missing.values * scale.factor + add.offset] <- NA

out.wx.data[1:rows,1:actual.columns,obs[i]] <- as.matrix(t.out)
}

## Update the status bar ##
if(!is.null(pb)){
cval <- pb$getVal()
Sys.sleep(0.000001)
setTkProgressBar(pb, cval+1, label=paste(round((cval+1)/increments*100, 0), "% done"))
}

		
loop.num <- loop.num + 1
unlink(c(out.temp,scale.offset.missingvals.temp))
}

#####################
## Print the units ##
if(return.units == TRUE){
	print(noquote(paste("Units of variable '", variable, "' are ", units, sep='')))
	}

###################
if(!is.null(pb)) { if(pb$getVal() == increments) {close(pb)} }
return(out.wx.data)
}

