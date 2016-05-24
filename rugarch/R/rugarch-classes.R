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

#----------------------------------------------------------------------------------
# Package Highest Level Virtual Class
#----------------------------------------------------------------------------------
setClass("rGARCH","VIRTUAL")
#----------------------------------------------------------------------------------
# univariate spec class
#----------------------------------------------------------------------------------
setClass("GARCHspec", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHspec",
		representation(model = "vector"),
		contains = "GARCHspec")
#----------------------------------------------------------------------------------
# univariate fit class
#----------------------------------------------------------------------------------
setClass("GARCHfit", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHfit",
		representation(fit = "vector",
				model = "vector"),
		contains = "GARCHfit")
#----------------------------------------------------------------------------------
# univariate filter class
#----------------------------------------------------------------------------------
setClass("GARCHfilter", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHfilter",
		representation(filter = "vector",
				model = "vector"),
		contains = "GARCHfilter")
#----------------------------------------------------------------------------------
# univariate forecast class (extends filter class)
#----------------------------------------------------------------------------------
setClass("GARCHforecast", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHforecast",
		representation(forecast = "vector",
				model = "vector"),
		contains = "GARCHforecast")
#----------------------------------------------------------------------------------
# univariate simulation class
#----------------------------------------------------------------------------------
setClass("GARCHsim", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHsim",
		representation(simulation = "vector",
				model = "vector",
				seed = "integer"),
		contains = "GARCHsim")
#----------------------------------------------------------------------------------
# univariate path simulation class
#----------------------------------------------------------------------------------
setClass("GARCHpath", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHpath",
		representation(path = "vector",
				model = "vector",
				seed = "integer"),
		contains = "GARCHpath")
#----------------------------------------------------------------------------------
# univariate garch roll class
#----------------------------------------------------------------------------------
setClass("GARCHroll", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHroll",
		representation(model = "vector",
				forecast = "vector"),
		contains = "GARCHroll")
#----------------------------------------------------------------------------------
# univariate garch parameter distribution (by simulation)
#----------------------------------------------------------------------------------
setClass("GARCHdistribution", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHdistribution",
		representation(dist = "vector",
				truecoef = "matrix",
				model = "vector"),
		contains =  "GARCHdistribution")
#----------------------------------------------------------------------------------
# univariate garch bootstrap forecast distribution
#----------------------------------------------------------------------------------
setClass("GARCHboot", contains = c("rGARCH", "VIRTUAL"))
setClass("uGARCHboot",
		representation(
				fseries = "matrix",
				fsigma = "matrix",
				bcoef = "data.frame",
				model = "vector",
				forc = "uGARCHforecast"),
		contains =  "GARCHboot")

#----------------------------------------------------------------------------------
# univariate garch test class
#----------------------------------------------------------------------------------
setClass("GARCHtests", contains = c("rGARCH", "VIRTUAL"))


#----------------------------------------------------------------------------------
# multiple spec/fit/filter/forecast garch methods (used in 2-stage extension multivariate models)
#----------------------------------------------------------------------------------
# Multiple Spec List Class
setClass("uGARCHmultispec", 
		representation(spec = "vector",
				type = "character"))

.validspeclist = function(object){
	all(unlist(lapply(object@spec, FUN = function(x) is(x, "uGARCHspec"))))
}

setValidity("uGARCHmultispec", .validspeclist)

# Multiple Fit ACD List Class
setClass("uGARCHmultifit", 
		representation(fit = "vector",
				desc = "vector"))

.validfitlist = function(object){
	all(unlist(lapply(object@fit, FUN = function(x) is(x, "uGARCHfit"))))
}

setValidity("uGARCHmultifit", .validfitlist)

# Multiple Fit ACD List Class
setClass("uGARCHmultifilter", 
		representation(filter = "vector",
				desc = "vector"))

.validfilterlist = function(object){
	all(unlist(lapply(object@filter, FUN = function(x) is(x, "uGARCHfilter"))))
}

setValidity("uGARCHmultifilter", .validfilterlist)

# Multiple Forecast ACD List Class
setClass("uGARCHmultiforecast", 
		representation(forecast = "vector",
				desc = "vector"))

.validforecastlist = function(object){
	all(unlist(lapply(object@forecast, FUN = function(x) is(x, "uGARCHforecast"))))
}

setValidity("uGARCHmultiforecast", .validforecastlist)


