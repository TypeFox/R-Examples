NCEP.gather.surface <-
function(variable, months.minmax, years.minmax, lat.minmax,
	lon.minmax,  reanalysis2=FALSE, return.units=TRUE, increments=NULL, pb=NULL) {


## Latitude and longitude should be given in decimal degrees ##
## 'months.minmax' must be given in the numeric format ##
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


##### Confirm that the variable selected exists at the level specified #####
if(any(c('air.sig995','lftx.sfc','lftx4.sfc','omega.sig995','pottmp.sig995','pr_wtr.eatm','pres.sfc','rhum.sig995','slp','mslp','uwnd.sig995','vwnd.sig995') == variable) == FALSE) { stop(paste("'",variable,"'"," not a valid variable name in reference to the surface.",sep='')) }

######################################################
## Specify which variables are not in Reanalysis II ##
not.in.reanalysis2 <- c('air.sig995','lftx.sfc','lftx4.sfc','omega.sig995','pottmp.sig995','rhum.sig995','slp','uwnd.sig995','vwnd.sig995')
not.in.reanalysis1 <- c('mslp')

#######################################################################################
#### Confirm that the variable selected exists in the Reanalysis dataset specified ####
if(reanalysis2 == TRUE && any(not.in.reanalysis2 == variable)) {
	reanalysis2 <- FALSE
	warning("This variable is not available at the surface in the Reanalysis II dataset.  Using Reanalysis I instead.")
}
if(reanalysis2 == FALSE && any(not.in.reanalysis1 == variable)) {
	reanalysis2 <- TRUE
	warning("This variable is not available at the surface in the Reanalysis I dataset.  Using Reanalysis II instead.")
}


#############################################
#############################################

## A lookup table for latitudes and longitudes
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

############################################################################################################################
## Figure out which latitudes and longitudes should be obtained such that the lat and lon input from the user is retained ##
## Latitudes ##
lat.minmax[1] <- floor(lat.minmax[1] / 2.5) * 2.5
lat.minmax[2] <- (ceiling(lat.minmax[2] / 2.5) * 2.5)
## Longitude ##
lon.minmax[1] <- floor(lon.minmax[1] / 2.5) * 2.5
lon.minmax[2] <- ifelse((ceiling(lon.minmax[2] / 2.5) * 2.5) > 357.5, 357.5, (ceiling(lon.minmax[2] / 2.5) * 2.5))



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
out.wx.data <- array(data=NA, dim=c(length(seq(lat.minmax[1],lat.minmax[2], by=2.5)), length(seq(lon.minmax[1],lon.minmax[2], by=2.5)), length(time.names)), 
	dimnames=list(rev(seq(lat.minmax[1],lat.minmax[2], by=2.5)),seq(lon.minmax[1],lon.minmax[2], by=2.5), time.names))

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
columns <- length(seq(lon.minmax[1],lon.minmax[2], by=2.5)) + 1
actual.columns <- columns-1
rows <- length(seq(lat.minmax[1],lat.minmax[2], by=2.5))


if(loop.num == 1){
### First query the data file online to determine the appropriate offset and scale factor to apply ##
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/surface/",variable,".",year,".nc.das", sep=''), scale.offset.missingvals.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/surface/",variable,".",year,".nc.das into a web browser to obtain an error message.", sep = ""))}
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
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
trying.out <- 1
fail <- 0
while(trying.out != 0){
trying.out <- try(download.file(paste("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/surface/",variable,".",year,".nc.ascii?",name,"[",beg.jdate,":",end.jdate,"][",lat.range[1],":",ifelse(length(lat.range) > 1, lat.range[2], lat.range[1]),"][",lon.range[1],":",ifelse(length(lon.range) > 1, lon.range[2], lon.range[1]),"]", sep=''), out.temp), silent=TRUE)
fail <- fail + 1
if(fail >= 5) {stop(paste("\nThere is a problem connecting to the NCEP database with the information provided.
	\nTry entering http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis",ifelse(reanalysis2 == TRUE, "2",""),"/surface/",variable,".",year,".nc.ascii?",name,"[",beg.jdate,":",end.jdate,"][",lat.range[1],":",ifelse(length(lat.range) > 1, lat.range[2], lat.range[1]),"][",lon.range[1],":",ifelse(length(lon.range) > 1, lon.range[2], lon.range[1]),"] into a web browser to obtain an error message.", sep=""))}
}


## Retrieve weather data from the temporary file ##
#### FOR VARIABLES AT OR NEAR THE SURFACE ####
outdata <- read.table(file=out.temp, sep=',', skip=12, header=FALSE, na.strings=missing.values, nrows=((end.jdate-beg.jdate)+1)*rows)


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

