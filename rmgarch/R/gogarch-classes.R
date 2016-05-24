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

# GO-GARCH Model (Factor/BEKK Family)
setClass("goGARCHspec",
		representation(model = "vector",
				umodel = "vector"),	
		contains = "mGARCHspec")

setClass("goGARCHfit",
		representation(mfit = "vector",
				model = "vector"),
		contains = "mGARCHfit")

setClass("goGARCHfilter",
		representation(mfilter = "vector",
				model = "vector"),
		contains = "mGARCHfilter")

setClass("goGARCHforecast",
		representation(mforecast = "vector",
				model = "vector"),
		contains = "mGARCHforecast")

setClass("goGARCHsim",
		representation(msim = "vector",
				model = "vector"),
		contains = "mGARCHsim")

setClass("goGARCHroll",
		representation(forecast = "vector",
				model = "vector"),
		contains = "mGARCHroll")

setClass("goGARCHfft",
		representation(dist = "vector", model = "vector"))
