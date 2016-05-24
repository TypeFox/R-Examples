#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


###################################################################################
#----------------------------------------------------------------------------------
# Multivariate Class Models
#----------------------------------------------------------------------------------
setClass("mGARCHspec", contains = c("GARCHspec"))
setClass("mGARCHfit", contains = c("GARCHfit"))
setClass("mGARCHfilter", contains = c("GARCHfilter"))
setClass("mGARCHforecast", contains = c("GARCHforecast"))
setClass("mGARCHsim", contains = c("GARCHsim"))
setClass("mGARCHroll", contains = c("GARCHroll"))
#----------------------------------------------------------------------------------

fScenario.default = function(){
	xmodel = list(asset.names = "", 
			assets = 0, 
			model = "dcc", 
			sim = 1, 
			roll = 0, 
			save.output = TRUE, 
			save.dir = "",
			save.name = "")
	return(list(scenario = list(), model = xmodel))
}

fMoments.default = function(){
	xmodel = list(asset.names = "", 
			assets = 0, 
			model = "dcc", 
			n.ahead = 1, 
			roll = 0, 
			save.output = TRUE, 
			save.dir = "",
			save.name = "")
	return(list(scenario = list(), model = xmodel))
}


setClass("fScenario", 
		representation(
				scenario = "vector", 
				model = "vector"),
		prototype = fScenario.default())

setClass("fMoments", 
		representation(
				moments = "vector", 
				model = "vector"),
		prototype = fMoments.default())