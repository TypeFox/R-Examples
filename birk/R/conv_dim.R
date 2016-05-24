# Created by Matthew A. Birk
# Converts between dimensions given a transition value
# Last updated: Jun 2015

#' Convert Dimensions of Measurement
#'
#' Converts between dimensions of measurement given a transition dimension (the dimension that "bridges" \code{x} and \code{y}, e.g. liters per second, lbs per acre). Note that 2 of the 3 measurements (\code{x}, \code{y}, or \code{trans}) must be defined to calculate the 3rd. See \code{\link{conv_unit_options}} for all options.
#'
#' This function supports all dimensions in \code{conv_unit_options} except for coordinates. The conversion values have been defined based primarily from international weight and measurement authorities (e.g. General Conference on Weights and Measures, International Committee for Weights and Measures, etc.). While much effort was made to make conversions as accurate as possible, you should check the accuracy of conversions to ensure that conversions are precise enough for your applications.
#'
#' @param x a numeric vector giving the measurement value in the first dimension.
#' @param x_unit the unit in which \code{x} was measured.
#' @param trans a numeric vector giving the measurement value in the transition dimension.
#' @param trans_unit the unit in which \code{trans} was measured.
#' @param y a numeric vector giving the measurement value in the second dimension.
#' @param y_unit the unit in which \code{y} was measured.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @note \describe{
#'	 \item{Duration}{Years are defined as 365.25 days and months are defined as 1/12 a year.}
#'	 \item{Energy}{cal is a thermochemical calorie (4.184 J) and Cal is 1000 cal (kcal or 4184 J).}
#'	 \item{Flow}{All gallon-based units are US gallons.}
#'	 \item{Mass}{All non-metric units are based on the avoirdupois system.}
#'	 \item{Power}{hp is mechanical horsepower, or 745.69 W.}
#'	 \item{Speed}{mach is calculated at sea level at 15 °C.}
#' }
#' @seealso \code{\link{conv_unit_options}}, \code{\link{conv_unit}}
#'
#' @examples
#' # How many minutes does it take to travel 100 meters at 3 feet per second?
#' conv_dim(x = 100, x_unit = "m", trans = 3, trans_unit = "ft_per_sec", y_unit = "min")
#' 
#' # How many degrees does the temperature increase with an increase in 4 kPa given 0.8 Celcius
#' # increase per psi?
#' conv_dim(x_unit = "C", trans = 0.8, trans_unit = "C_per_psi", y = 4, y_unit = "kPa")
#' 
#' # Find the densities given volume and mass measurements.
#' conv_dim(x = c(60, 80), x_unit = "ft3", trans_unit = "kg_per_l", y = c(6e6, 4e6), y_unit = "g")
#' 
#' @encoding UTF-8
#' @export

conv_dim = function(x, x_unit, trans, trans_unit, y, y_unit){
	
	if(trans_unit == 'grav'){
		orig_trans_unit = trans_unit
		trans_unit = 'm_per_sec2'
		trans = conv_unit(trans,orig_trans_unit,trans_unit)
	}
	if(trans_unit == 'Sv'){
		orig_trans_unit = trans_unit
		trans_unit = 'l_per_sec'
		trans = conv_unit(trans,orig_trans_unit,trans_unit)
	}
	if(length(grep('_per_sec2', trans_unit)) == 1) trans_unit = gsub('_per_sec2', '.per.sec_per_sec', trans_unit)
	if(trans_unit %in% c('uW', 'mW', 'W', 'kW', 'MW', 'GW', 'hp')){
		orig_trans_unit = trans_unit
		trans_unit = 'erg_per_sec'
		trans = conv_unit(trans, orig_trans_unit, trans_unit)
	}
	if(trans_unit %in% c('kph', 'mph', 'knot', 'mach', 'light')){
		orig_trans_unit = trans_unit
		trans_unit = 'm_per_sec'
		trans = conv_unit(trans, orig_trans_unit, trans_unit)
	}
	
	units = unlist(strsplit(trans_unit, split = '_per_'))
	units = gsub('.per.', '_per_', units) # for acceleration only
	dims = sapply(units, function(i) names(conv_unit_options)[sapply(conv_unit_options, function(j) i %in% j)])
	if(any(sapply(dims, length) == 0)) stop('the \'trans_unit\' argument is not an acceptable unit. It must either be listed in conv_unit_options or a combination of units from conv_unit_options separated by \"_per_\"')
	x_dim = names(conv_unit_options)[sapply(conv_unit_options, function(j) x_unit %in% j)]
	if(length(x_dim) == 0) stop('the \'x_unit\' argument is not an acceptable unit')
	y_dim = names(conv_unit_options)[sapply(conv_unit_options, function(j) y_unit %in% j)]
	if(length(y_dim) == 0) stop('the \'y_unit\' argument is not an acceptable unit')
	if(!missing('x')) x = conv_unit(x, x_unit, names(which(dims == x_dim)))
	if(!missing('y')) y = conv_unit(y, y_unit, names(which(dims == y_dim)))
	
if(missing(y)){
	if(which(dims == x_dim) == 1){
		y = x/trans
		y = conv_unit(y, names(dims[which(dims == y_dim)]), y_unit)
	}else{
		y = x*trans
		y = conv_unit(y, names(dims[which(dims == y_dim)]), y_unit)
	}
	return(y)
}
	
if(missing(trans)){
	if(all(c(x_dim, y_dim) == dims)) trans = x/y else trans = y/x
	if(exists('orig_trans_unit', inherits = F)) trans = conv_unit(trans, trans_unit, orig_trans_unit)
	return(trans)
}
	
if(missing(x)){
	if(which(dims == y_dim) == 1){
		x = y/trans
		x = conv_unit(x, names(dims[which(dims == x_dim)]), x_unit)
	}else{
		x = y*trans
		x = conv_unit(x, names(dims[which(dims == x_dim)]), x_unit)
	}
	return(x)
}
}