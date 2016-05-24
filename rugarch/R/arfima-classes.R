#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


setClass("ARFIMA", contains = c("rGARCH", "VIRTUAL"))

# ARFIMA Model
setClass("ARFIMAspec",
		representation(model = "vector"),
		contains = "ARFIMA")

setClass("ARFIMAfilter",
		representation(filter = "vector",
				model = "vector"),
		contains = "ARFIMA")


setClass("ARFIMAfit",
		representation(fit = "vector",
				model = "vector"),
		contains = "ARFIMA")

setClass("ARFIMAsim",
		representation(simulation = "vector",
				model = "vector",
				seed = "integer"),
		contains = "ARFIMA")

setClass("ARFIMAforecast",
		representation(forecast = "vector",
				model = "vector"),
		contains = "ARFIMA")

setClass("ARFIMApath",
		representation(path = "vector",
				model = "vector",
				seed = "integer"),
		contains = "ARFIMA")

setClass("ARFIMAroll",
		representation(model = "vector",
				forecast = "vector"),
		contains = "ARFIMA")


setClass("ARFIMAdistribution",
		representation(dist = "vector",
				truecoef = "matrix",
				model = "vector"),
		contains =  "ARFIMA")

# Have not yet implemented the ARFIMA bootstrap and not likely to in near future
setClass("ARFIMAboot",
		representation(
				fseries = "matrix",
				bcoef = "data.frame",
				model = "vector",
				forc = "ARFIMAforecast"),
		contains =  "ARFIMA")
#----------------------------------------------------------------------------------
# multiple spec/fit/filter/forecast ARFIMA methods
#----------------------------------------------------------------------------------

# Multiple Spec List Class
setClass("ARFIMAmultispec", 
		representation(spec = "vector",
				type = "character"),
		contains = "ARFIMA")

.validarfimaspeclist = function(object){
	all(unlist(lapply(object@spec, FUN = function(x) is(x, "ARFIMAspec"))))
}

setValidity(Class = "ARFIMAmultispec", method = .validarfimaspeclist)

# Multiple Fit ACD List Class
setClass("ARFIMAmultifit", 
		representation(fit = "vector",
				desc = "vector"),
		contains = "ARFIMA")

.validarfimafitlist = function(object){
	all(unlist(lapply(object@fit, FUN = function(x) is(x, "ARFIMAfit"))))
}

setValidity("ARFIMAmultifit", .validarfimafitlist)

# Multiple Fit ACD List Class
setClass("ARFIMAmultifilter", 
		representation(filter = "vector",
				desc = "vector"),
		contains = "ARFIMA")

.validarfimafilterlist = function(object){
	all(unlist(lapply(object@filter, FUN = function(x) is(x, "ARFIMAfilter"))))
}

setValidity("ARFIMAmultifilter", .validarfimafilterlist)


# Multiple Forecast ACD List Class
setClass("ARFIMAmultiforecast", 
		representation(forecast = "vector",
				desc = "vector"),
		contains = "ARFIMA")

.validarfimaforecastlist = function(object){
	all(unlist(lapply(object@forecast, FUN = function(x) is(x, "ARFIMAforecast"))))
}

setValidity("ARFIMAmultiforecast", .validarfimaforecastlist)