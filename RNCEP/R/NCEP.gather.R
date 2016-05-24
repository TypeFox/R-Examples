NCEP.gather <-
function(variable, level, months.minmax, years.minmax, 
	lat.southnorth, lon.westeast, reanalysis2=FALSE, return.units=TRUE, status.bar=TRUE){

## Make sure that the reference system has been specified ##
if(is.null(level)) { stop("One of 'surface', 'gaussian', or a numeric pressure level must be given for 'level'") }

## Make sure that the user only specified one reference system ##
if(length(level) > 1) { stop("Cannot access multiple reference systems in a single function call") }

## Make sure that the reference has been specified properly ##
if(is.numeric(level) == FALSE) { 
	if(level %in% c('surface','gaussian') == FALSE) { stop("level must be one of 'gaussian', 'surface' or a numeric pressure level") }
	}

## Make sure that the data exist for the dates given ##
if(reanalysis2 == TRUE && years.minmax[1] < 1979) {stop("The datetimes specified are out of range for the Reanalysis 2 dataset.")}
if(years.minmax[1] < 1948) {stop("The datetimes specified are out of range.")}

## Convert negative longitude values to positive ones ##
lon.westeast[1] <- ifelse(lon.westeast[1] < 0, 360+lon.westeast[1], lon.westeast[1])
lon.westeast[length(lon.westeast)] <- ifelse(lon.westeast[length(lon.westeast)] < 0, 360+lon.westeast[length(lon.westeast)], lon.westeast[length(lon.westeast)])


## Determine whether data will need to be obtained from both sides of the Prime Meridian ##
if(lon.westeast[1] > lon.westeast[length(lon.westeast)]) { cross.prime <- TRUE } else { cross.prime <- FALSE }

## If a status bar is desired, describe the status bar parameters ##
if(status.bar){
	#importFrom(tcltk,tkProgressBar)
	#require(tcltk)
if(cross.prime){
	tlength <- length(seq(years.minmax[1],years.minmax[length(years.minmax)]))*2
	pb <- tkProgressBar(title="Total progress", min = 0,
		max=tlength, width=300) } else {
	tlength <- length(seq(years.minmax[1],years.minmax[length(years.minmax)]))
	pb <- tkProgressBar(title="Total progress", min = 0,
		max=tlength, width=300) }
	} else { 
		tlength <- NULL
		pb <- NULL 
		}


## Apply the appropriate function for queries that DO NOT cross the Prime Meridian ##
if(cross.prime == FALSE){
	if(is.numeric(level)){
		out <- NCEP.gather.pressure(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = lon.westeast, pressure = level, reanalysis2 = reanalysis2, return.units = return.units,
			pb=pb, increments=tlength)
	} else
	if(level == 'surface'){
		out <- NCEP.gather.surface(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = lon.westeast, reanalysis2 = reanalysis2, return.units = return.units, 
			pb=pb, increments=tlength)
	} else
	if(level == 'gaussian'){
		out <- NCEP.gather.gaussian(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = lon.westeast, reanalysis2 = reanalysis2, return.units = return.units,
			pb=pb, increments=tlength)
	}
} else

## Apply the appropriate function for queries that DO cross the Prime Meridian ##
if(cross.prime == TRUE){
	if(is.numeric(level)){
		out.west <- NCEP.gather.pressure(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(lon.westeast[1], 357.5), pressure = level, reanalysis2 = reanalysis2, 
			return.units = return.units, pb=pb, increments=tlength)
		out.east <- NCEP.gather.pressure(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(0,lon.westeast[2]), pressure = level, reanalysis2 = reanalysis2, 
			return.units = return.units, pb=pb, increments=tlength)
		out <- NCEP.bind(data.west = out.west, data.east = out.east)
	} else
	if(level == 'surface'){
		out.west <- NCEP.gather.surface(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(lon.westeast[1], 357.5), reanalysis2 = reanalysis2, 
			return.units = return.units, pb=pb, increments=tlength)
		out.east <- NCEP.gather.surface(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(0,lon.westeast[2]), reanalysis2 = reanalysis2,
			return.units = return.units, pb=pb, increments=tlength)
		out <- NCEP.bind(data.west = out.west, data.east = out.east)
	} else
	if(level == 'gaussian'){
		out.west <- NCEP.gather.gaussian(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(lon.westeast[1],358.125), reanalysis2 = reanalysis2, 
			return.units = return.units, pb=pb, increments=tlength)
		out.east <- NCEP.gather.gaussian(variable=variable, months.minmax = months.minmax, years.minmax = years.minmax, 
			lat.minmax = lat.southnorth, lon.minmax = c(0,lon.westeast[2]), reanalysis2 = reanalysis2, 
			return.units = return.units, pb=pb, increments=tlength)
		out <- NCEP.bind(data.west = out.west, data.east = out.east)
	}
}

return(out)

}

