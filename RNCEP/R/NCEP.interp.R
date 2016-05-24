NCEP.interp <-
function(variable, level, lat, lon, dt,	reanalysis2=FALSE, 
	interpolate.space=TRUE, interpolate.time=TRUE, keep.unpacking.info=FALSE,
	return.units=TRUE, interp='linear', p=1, status.bar=TRUE){

## Make sure that the reference system has been specified ##
if(is.null(level)) { stop("One of 'surface', 'gaussian', or a numeric pressure level must be given for 'level'") }

## Make sure that the user only specified one reference system ##
if(length(level) > 1 && is.numeric(level) == FALSE) { stop("Cannot access multiple levels (except pressure levels) in a single function call") }

## Make sure that the reference has been specified properly ##
if(is.numeric(level) == FALSE) { 
	if(level %in% c('surface','gaussian') == FALSE) { stop("level must be one of 'gaussian', 'surface' or a numeric pressure level") }
	}

## Make sure that the data exist for the dates given ##
if(reanalysis2 == TRUE && as.numeric(format(strptime(dt, "%Y-%m-%d %H:%M:%S"), "%Y")) < 1979) {stop("The datetime specified is out of range for the Reanalysis 2 dataset.")}
if(any(as.numeric(format(strptime(dt, "%Y-%m-%d %H:%M:%S"), "%Y")) < 1948)) {stop("The datetime specified is out of range.")}

## Make sure that keep.unpacking.info == FALSE when it should be ##
if(length(unique(variable)) > 1 && keep.unpacking.info==TRUE) {
	keep.unpacking.info=FALSE
	warning("keep.unpacking.info was changed to FALSE because multiple variables were queried")
	}
if(length(unique(reanalysis2)) > 1 && keep.unpacking.info==TRUE) {
	keep.unpacking.info=FALSE
	warning("keep.unpacking.info was changed to FALSE because data were obtained from both Reanalysis I and Reanalysis II datasets")
	}


## Make sure that the latitudes and longitudes are within the valid range ##
## Create 'is.between' function ##
is.between <- function(x, a, b) {
(x - a)  *  (b - x) >= 0
}
if(any(is.between(lat, -90, 90) == FALSE)) { stop("Latitudes must be between -90 and 90") }
if(any(is.between(lon, -360, 360) == FALSE)) { stop("Longitudes must be between -360 and 360") }

## Apply the correct function ##
if(is.numeric(level)){
	out <- NCEP.interp.pressure(variable=variable, lat=lat, lon=lon, dt=dt, pressure=level, reanalysis2=reanalysis2, interpolate.space=interpolate.space, interpolate.time=interpolate.time, keep.unpacking.info=keep.unpacking.info, return.units=return.units, interp=interp, p=p, status.bar=status.bar)
} else
if(level == 'surface'){
	out <- NCEP.interp.surface(variable=variable, lat=lat, lon=lon, dt=dt, reanalysis2=reanalysis2, interpolate.space=interpolate.space, interpolate.time=interpolate.time, keep.unpacking.info=keep.unpacking.info, return.units=return.units, interp=interp, p=p, status.bar=status.bar)
} else
if(level == 'gaussian'){
	out <- NCEP.interp.gaussian(variable=variable, lat=lat, lon=lon, dt=dt, reanalysis2=reanalysis2, interpolate.space=interpolate.space, interpolate.time=interpolate.time, keep.unpacking.info=keep.unpacking.info, return.units=return.units, interp=interp, p=p, status.bar=status.bar)
}

return(out)

}

