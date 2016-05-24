# Created by Matthew A. Birk
# Converts common units for a variety of dimensions
# Last updated: Feb 2016

.conversions=data.frame(dim = character(0), unit = character(0), std = numeric(0))
.conversions = rbind(.conversions,
	
	data.frame(dim = 'acceleration', unit = 'mm_per_sec2', std = 100),
	data.frame(dim = 'acceleration', unit = 'cm_per_sec2', std = 10),
	data.frame(dim = 'acceleration', unit = 'm_per_sec2', std = 1),
	data.frame(dim = 'acceleration', unit = 'km_per_sec2', std = 1e-3),
	data.frame(dim = 'acceleration', unit = 'grav', std = 1/9.80665),
	data.frame(dim = 'acceleration', unit = 'inch_per_sec2', std = 100/2.54),
	data.frame(dim = 'acceleration', unit = 'ft_per_sec2', std = 100/2.54/12),
	data.frame(dim = 'acceleration', unit = 'mi_per_sec2', std = 100/2.54/12/5280),
	data.frame(dim = 'acceleration', unit = 'kph_per_sec', std = 1e-3*3600),
	data.frame(dim = 'acceleration', unit = 'mph_per_sec', std = 100/2.54/12/5280*3600),
	
	data.frame(dim = 'angle', unit = 'degree', std = 360),
	data.frame(dim = 'angle', unit = 'radian', std = 2 * pi),
	data.frame(dim = 'angle', unit = 'grad', std = 400),
	data.frame(dim = 'angle', unit = 'arcmin', std = 21600),
	data.frame(dim = 'angle', unit = 'arcsec', std = 1296000),
	data.frame(dim = 'angle', unit = 'turn', std = 1),

	data.frame(dim = 'area', unit = 'nm2', std = 1e18),
	data.frame(dim = 'area', unit = 'um2', std = 1e12),
	data.frame(dim = 'area', unit = 'mm2', std = 1e6),
	data.frame(dim = 'area', unit = 'cm2', std = 1e4),
	data.frame(dim = 'area', unit = 'm2', std = 1),
	data.frame(dim = 'area', unit = 'hectare', std = 1e-4),
	data.frame(dim = 'area', unit = 'km2', std = 1e-6),
	data.frame(dim = 'area', unit = 'inch2', std = (100/2.54)^2),
	data.frame(dim = 'area', unit = 'ft2', std = (100/2.54/12)^2),
	data.frame(dim = 'area', unit = 'yd2', std = (100/2.54/36)^2),
	data.frame(dim = 'area', unit = 'acre', std = (100/2.54/12)^2/43560),
	data.frame(dim = 'area', unit = 'mi2', std = (100/2.54/(12*5280))^2),
	data.frame(dim = 'area', unit = 'naut_mi2', std = 1/3429904),

	data.frame(dim = 'coordinate', unit = 'dec_deg', std = NA),
	data.frame(dim = 'coordinate', unit = 'deg_dec_min', std = NA),
	data.frame(dim = 'coordinate', unit = 'deg_min_sec', std = NA),

	data.frame(dim = 'count', unit = 'nmol', std = 1e9),
	data.frame(dim = 'count', unit = 'umol', std = 1e6),
	data.frame(dim = 'count', unit = 'mmol', std = 1e3),
	data.frame(dim = 'count', unit = 'mol', std = 1),

	data.frame(dim = 'duration', unit = 'nsec', std = 1e9),
	data.frame(dim = 'duration', unit = 'usec', std = 1e6),
	data.frame(dim = 'duration', unit = 'msec', std = 1e3),
	data.frame(dim = 'duration', unit = 'sec', std = 1),
	data.frame(dim = 'duration', unit = 'min', std = 1/60),
	data.frame(dim = 'duration', unit = 'hr', std = 1/3600),
	data.frame(dim = 'duration', unit = 'day', std = 1/86400),
	data.frame(dim = 'duration', unit = 'wk', std = 1/604800),
	data.frame(dim = 'duration', unit = 'mon', std = 1/(86400*365.25/12)),
	data.frame(dim = 'duration', unit = 'yr', std = 1/(86400*365.25)),
	data.frame(dim = 'duration', unit = 'dec', std = 1/(86400*3652.5)),
	data.frame(dim = 'duration', unit = 'cen', std = 1/(86400*36525)),
	data.frame(dim = 'duration', unit = 'mil', std = 1/(86400*365250)),
	data.frame(dim = 'duration', unit = 'Ma', std = 1/(86400*365250000)),

	data.frame(dim = 'energy', unit = 'J', std = 1),
	data.frame(dim = 'energy', unit = 'kJ', std = 1e-3),
	data.frame(dim = 'energy', unit = 'erg', std = 1e7),
	data.frame(dim = 'energy', unit = 'cal', std = 1/4.184),
	data.frame(dim = 'energy', unit = 'Cal', std = 1/4184),
	data.frame(dim = 'energy', unit = 'Wsec', std = 1),
	data.frame(dim = 'energy', unit = 'kWh', std = 1/3.6e6),
	data.frame(dim = 'energy', unit = 'MWh', std = 1/3.6e9),
	data.frame(dim = 'energy', unit = 'BTU', std = 1/1055.05585262),

	data.frame(dim = 'flow', unit = 'ml_per_sec', std = 1e3),
	data.frame(dim = 'flow', unit = 'ml_per_min', std = 1e3*60),
	data.frame(dim = 'flow', unit = 'ml_per_hr', std = 1e3*3600),
	data.frame(dim = 'flow', unit = 'l_per_sec', std = 1),
	data.frame(dim = 'flow', unit = 'l_per_min', std = 60),
	data.frame(dim = 'flow', unit = 'l_per_hr', std = 3600),
	data.frame(dim = 'flow', unit = 'm3_per_sec', std = 1e-3),
	data.frame(dim = 'flow', unit = 'm3_per_min', std = 1e-3*60),
	data.frame(dim = 'flow', unit = 'm3_per_hr', std = 1e-3*3600),
	data.frame(dim = 'flow', unit = 'Sv', std = 1e-9),
	data.frame(dim = 'flow', unit = 'gal_per_sec', std = 1/3.785411784),
	data.frame(dim = 'flow', unit = 'gal_per_min', std = 1/3.785411784*60),
	data.frame(dim = 'flow', unit = 'gal_per_hr', std = 1/3.785411784*3600),
	data.frame(dim = 'flow', unit = 'ft3_per_sec', std = 61.0237440947323/1728),
	data.frame(dim = 'flow', unit = 'ft3_per_min', std = 61.0237440947323/1728*60),
	data.frame(dim = 'flow', unit = 'ft3_per_hr', std = 61.0237440947323/1728*3600),
	
	data.frame(dim = 'length', unit = 'angstrom', std = 1e10),
	data.frame(dim = 'length', unit = 'nm', std = 1e9),
	data.frame(dim = 'length', unit = 'um', std = 1e6),
	data.frame(dim = 'length', unit = 'mm', std = 1e3),
	data.frame(dim = 'length', unit = 'cm', std = 100),
	data.frame(dim = 'length', unit = 'dm', std = 10),
	data.frame(dim = 'length', unit = 'm', std = 1),
	data.frame(dim = 'length', unit = 'km', std = 1e-3),
	data.frame(dim = 'length', unit = 'inch', std = 100/2.54),
	data.frame(dim = 'length', unit = 'ft', std = 100/2.54/12),
	data.frame(dim = 'length', unit = 'yd', std = 100/2.54/36),
	data.frame(dim = 'length', unit = 'fathom', std = 100/2.54/72),
	data.frame(dim = 'length', unit = 'mi', std = 100/2.54/12/5280),
	data.frame(dim = 'length', unit = 'naut_mi', std = 1/1852),
	data.frame(dim = 'length', unit = 'au', std = 1/149597870700),
	data.frame(dim = 'length', unit = 'light_yr', std = 1/9460730472580800),
	data.frame(dim = 'length', unit = 'parsec', std = 1/149597870700/(6.48e5/pi)),
	data.frame(dim = 'length', unit = 'point', std = 100/2.54*72),
	
	data.frame(dim = 'mass', unit = 'ug', std = 1e6),
	data.frame(dim = 'mass', unit = 'mg', std = 1e3),
	data.frame(dim = 'mass', unit = 'g', std = 1),
	data.frame(dim = 'mass', unit = 'kg', std = 1e-3),
	data.frame(dim = 'mass', unit = 'Pg', std = 1e-15),
	data.frame(dim = 'mass', unit = 'carat', std = 5),
	data.frame(dim = 'mass', unit = 'metric_ton', std = 1e-6),
	data.frame(dim = 'mass', unit = 'oz', std = 1/28.349523125),
	data.frame(dim = 'mass', unit = 'lbs', std = 2.20462234e-3),
	data.frame(dim = 'mass', unit = 'short_ton', std = 1/907184.74),
	data.frame(dim = 'mass', unit = 'long_ton', std = 1/1.016e6),
	data.frame(dim = 'mass', unit = 'stone', std = 2.20462234e-3/14),
	
	data.frame(dim = 'power', unit = 'uW', std = 1e6),
	data.frame(dim = 'power', unit = 'mW', std = 1e3),
	data.frame(dim = 'power', unit = 'W', std = 1),
	data.frame(dim = 'power', unit = 'kW', std = 1e-3),
	data.frame(dim = 'power', unit = 'MW', std = 1e-6),
	data.frame(dim = 'power', unit = 'GW', std = 1e-9),
	data.frame(dim = 'power', unit = 'erg_per_sec', std = 1e7),
	data.frame(dim = 'power', unit = 'cal_per_sec', std = 1/4.184),
	data.frame(dim = 'power', unit = 'cal_per_hr', std = 1/4.184*3600),
	data.frame(dim = 'power', unit = 'Cal_per_sec', std = 1/4184),
	data.frame(dim = 'power', unit = 'Cal_per_hr', std = 1/4184*3600),
	data.frame(dim = 'power', unit = 'BTU_per_sec', std = 1/1055.05585262),
	data.frame(dim = 'power', unit = 'BTU_per_hr', std = 1/1055.05585262*3600),
	data.frame(dim = 'power', unit = 'hp', std = 1/745.69),
	
	data.frame(dim = 'pressure', unit = 'uatm', std = 1e6),
	data.frame(dim = 'pressure', unit = 'atm', std = 1),
	data.frame(dim = 'pressure', unit = 'Pa', std = 101325),
	data.frame(dim = 'pressure', unit = 'hPa', std = 1013.25),
	data.frame(dim = 'pressure', unit = 'kPa', std = 101.325),
	data.frame(dim = 'pressure', unit = 'torr', std = 760),
	data.frame(dim = 'pressure', unit = 'mmHg', std = 760),
	data.frame(dim = 'pressure', unit = 'inHg', std = 1/(3386.389/101325)),
	data.frame(dim = 'pressure', unit = 'mbar', std = 1013.25),
	data.frame(dim = 'pressure', unit = 'bar', std = 1.01325),
	data.frame(dim = 'pressure', unit = 'dbar', std = 10.1325),
	data.frame(dim = 'pressure', unit = 'psi', std = 14.69594877551),
	
	data.frame(dim = 'speed', unit = 'mm_per_sec', std = 1e3),
	data.frame(dim = 'speed', unit = 'cm_per_sec', std = 100),
	data.frame(dim = 'speed', unit = 'm_per_sec', std = 1),
	data.frame(dim = 'speed', unit = 'km_per_sec', std = 1e-3),
	data.frame(dim = 'speed', unit = 'inch_per_sec', std = 100/2.54),
	data.frame(dim = 'speed', unit = 'ft_per_sec', std = 100/2.54/12),
	data.frame(dim = 'speed', unit = 'kph', std = 1e-3*3600),
	data.frame(dim = 'speed', unit = 'mph', std = 100/2.54/12/5280*3600),
	data.frame(dim = 'speed', unit = 'km_per_day', std = 1e-3*3600*24),
	data.frame(dim = 'speed', unit = 'mi_per_day', std = 100/2.54/12/5280*3600*24),
	data.frame(dim = 'speed', unit = 'knot', std = 1/1852*3600),
	data.frame(dim = 'speed', unit = 'mach', std = 1/340.3),
	data.frame(dim = 'speed', unit = 'light', std = 1/299792458),
	
	data.frame(dim = 'temperature', unit = 'C', std = NA),
	data.frame(dim = 'temperature', unit = 'F', std = NA),
	data.frame(dim = 'temperature', unit = 'K', std = NA),
	data.frame(dim = 'temperature', unit = 'R', std = NA),
	
	data.frame(dim = 'volume', unit = 'ul', std = 1e6),
	data.frame(dim = 'volume', unit = 'ml', std = 1e3),
	data.frame(dim = 'volume', unit = 'dl', std = 10),
	data.frame(dim = 'volume', unit = 'l', std = 1),
	data.frame(dim = 'volume', unit = 'cm3', std = 1e3),
	data.frame(dim = 'volume', unit = 'dm3', std = 1),
	data.frame(dim = 'volume', unit = 'm3', std = 1e-3),
	data.frame(dim = 'volume', unit = 'km3', std = 1e-12),
	data.frame(dim = 'volume', unit = 'us_tsp', std = 1/3.785411784*768),
	data.frame(dim = 'volume', unit = 'us_tbsp', std = 1/3.785411784*256),
	data.frame(dim = 'volume', unit = 'us_oz', std = 1/3.785411784*128),
	data.frame(dim = 'volume', unit = 'us_cup', std = 1/3.785411784*16),
	data.frame(dim = 'volume', unit = 'us_pint', std = 1/3.785411784*8),
	data.frame(dim = 'volume', unit = 'us_quart', std = 1/3.785411784*4),
	data.frame(dim = 'volume', unit = 'us_gal', std = 1/3.785411784),
	data.frame(dim = 'volume', unit = 'inch3', std = 61.0237440947323),
	data.frame(dim = 'volume', unit = 'ft3', std = 61.0237440947323/12^3),
	data.frame(dim = 'volume', unit = 'mi3', std = 61.0237440947323/(12*5280)^3),
	data.frame(dim = 'volume', unit = 'imp_tsp', std = 1/4.54609*768),
	data.frame(dim = 'volume', unit = 'imp_tbsp', std = 1/4.54609*256),
	data.frame(dim = 'volume', unit = 'imp_oz', std = 1/4.54609*160),
	data.frame(dim = 'volume', unit = 'imp_cup', std = 1/4.54609*16),
	data.frame(dim = 'volume', unit = 'imp_pint', std = 1/4.54609*8),
	data.frame(dim = 'volume', unit = 'imp_quart', std = 1/4.54609*4),
	data.frame(dim = 'volume', unit = 'imp_gal', std = 1/4.54609)
)

conv_unit_options = lapply(split(.conversions$unit, .conversions$dim, drop = T), as.character)


#' Convert Units of Measurement
#'
#' Converts common units of measurement for a variety of dimensions. See \code{\link{conv_unit_options}} for all options.
#'
#' \describe{
#'	 \item{Acceleration}{mm_per_sec2, cm_per_sec2, m_per_sec2, km_per_sec2, grav, inch_per_sec2, ft_per_sec2, mi_per_sec2, kph_per_sec, mph_per_sec}
#'	 \item{Angle}{degree, radian, grad, arcmin, arcsec, turn}
#'	 \item{Area}{nm2, um2, mm2, cm2, m2, hectare, km2, inch2, ft2, yd2, acre, mi2, naut_mi2}
#'	 \item{Coordinate}{dec_deg, deg_dec_min, deg_min_sec (see note)}
#'	 \item{Count}{nmol, umol, mmol, mol}
#'	 \item{Duration}{nsec, usec, msec, sec, min, hr, day, wk, mon, yr, dec, cen, mil, Ma}
#'	 \item{Energy}{J, kJ, erg, cal, Cal, Wsec, kWh, MWh, BTU}
#'	 \item{Flow}{ml_per_sec, ml_per_min, ml_per_hr, l_per_sec, l_per_min, l_per_hr, m3_per_sec, m3_per_min, m3_per_hr, gal_per_sec, gal_per_min, gal_per_hr, ft3_per_sec, ft3_per_min, ft3_per_hr, Sv}
#'	 \item{Length}{angstrom, nm, um, mm, cm, dm, m, km, inch, ft, yd, fathom, mi, naut_mi, au, light_yr, parsec, point}
#'	 \item{Mass}{ug, mg, g, kg, Pg, carat, metric_ton, oz, lbs, short_ton, long_ton, stone}
#'	 \item{Power}{uW, mW, W, kW, MW, GW, erg_per_sec, cal_per_sec, cal_per_hr, Cal_per_sec, Cal_per_hr, BTU_per_sec, BTU_per_hr, hp}
#'	 \item{Pressure}{uatm, atm, Pa, hPa, kPa, torr, mmHg, inHg, mbar, bar, dbar, psi}
#'	 \item{Speed}{mm_per_sec, cm_per_sec, m_per_sec, km_per_sec, inch_per_sec, ft_per_sec, kph, mph, km_per_day, mi_per_day, knot, mach, light}
#'	 \item{Temperature}{C, F, K, R}
#'	 \item{Volume}{ul, ml, dl, l, cm3, dm3, m3, km3, us_tsp, us_tbsp, us_oz, us_cup, us_pint, us_quart, us_gal, inch3, ft3, mi3, imp_tsp, imp_tbsp, imp_oz, imp_cup, imp_pint, imp_quart, imp_gal}
#' }
#' The conversion values have been defined based primarily from international weight and measurement authorities (e.g. General Conference on Weights and Measures, International Committee for Weights and Measures, etc.). While much effort was made to make conversions as accurate as possible, you should check the accuracy of conversions to ensure that conversions are precise enough for your applications.
#'
#' @param x a numeric vector giving the measurement value in its original units.
#' @param from the unit in which the measurement was made.
#' @param to the unit to which the measurement is to be converted.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @note \describe{
#'	 \item{Duration}{Years are defined as 365.25 days and months are defined as 1/12 a year.}
#'	 \item{Coordinate}{Values must be entered as a string with one space between subunits (e.g. 70° 33’ 11” = "70 33 11").}
#'	 \item{Energy}{cal is a thermochemical calorie (4.184 J) and Cal is 1000 cal (kcal or 4184 J).}
#'	 \item{Flow}{All gallon-based units are US gallons.}
#'	 \item{Mass}{All non-metric units are based on the avoirdupois system.}
#'	 \item{Power}{hp is mechanical horsepower, or 745.69 W.}
#'	 \item{Speed}{mach is calculated at sea level at 15 °C.}
#' }
#' @seealso \code{\link{conv_unit_options}}, \code{\link{conv_dim}}
#'
#' @examples
#' conv_unit(2.54, "cm", "inch") # Result = 1 inch
#' 
#' conv_unit(seq(1, 10), "kg", "short_ton") # A vector of measurement values can be converted
#' 
#' # Convert 1, 10, and 100 meters to all other length units
#' sapply(conv_unit_options$length, function(x) conv_unit(c(1, 10, 100), "m", x))
#' 
#' conv_unit("33 1 1", "deg_min_sec", "dec_deg")
#' 
#' conv_unit(c("101 44.32","3 19.453"), "deg_dec_min", "deg_min_sec")
#'
#' @encoding UTF-8
#' @export

conv_unit=function(x,from,to)
{
	unit=std=NULL
	if(nrow(subset(.conversions,unit==from,dim))==0) stop('the \'from\' argument is not an acceptable unit.')
	if(nrow(subset(.conversions,unit==to,dim))==0) stop('the \'to\' argument is not an acceptable unit.')
	if(subset(.conversions,unit==from,dim)!=subset(.conversions,unit==to,dim)) stop('these units cannot be converted because they are of different dimensions. Try using conv_dim().')
	if((from=='C' | from=='F' | from=='K' | from=='R') & (to=='C' | to=='F' | to=='K' | to=='R'))
	{
		frzC=0.01
		frzF=32.018
		frzK=273.16
		frzR=491.688
		boilC=99.9839
		boilF=211.97102
		boilK=373.1339
		boilR=671.64102
		rangeC=boilC-frzC
		rangeF=boilF-frzF
		rangeK=boilK-frzK
		rangeR=boilR-frzR
		prop=(x-get(paste('frz',from,sep='')))/get(paste('range',from,sep=''))
		return(prop*get(paste('range',to,sep=''))+get(paste('frz',to,sep='')))
	}
	if((from=='dec_deg' | from=='deg_dec_min' | from=='deg_min_sec') & (to=='dec_deg' | to=='deg_dec_min' | to=='deg_min_sec'))
	{
		if(from=='dec_deg') secs=as.numeric(x)*3600
		if(from=='deg_dec_min') secs=lapply(split(as.numeric(unlist(strsplit(x,' ')))*c(3600,60),f=rep(1:length(x),each=2)),sum)
		if(from=='deg_min_sec') secs=lapply(split(as.numeric(unlist(strsplit(x,' ')))*c(3600,60,1),f=rep(1:length(x),each=3)),sum)
		if(to=='dec_deg') return(as.character(lapply(secs,function(y) y/3600)))
		if(to=='deg_dec_min') return(paste(lapply(secs,function(y) y%/%3600),lapply(secs,function(y) y%%3600/60)))
		if(to=='deg_min_sec') return(paste(lapply(secs,function(y) y%/%3600),lapply(secs,function(y) y%%3600%/%60),lapply(secs,function(y) y%%3600%%60)))
	}
	value=x/subset(.conversions,unit==from,std,drop=T)
	return(value*subset(.conversions,unit==to,std,drop=T))
}