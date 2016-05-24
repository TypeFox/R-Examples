# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MareyMap.
#
# MareyMap is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MareyMap is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MareyMap; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


#--- accessors
setMethod("rates", "Interpolation", function(object) object@rates)
setMethod("visible", "Interpolation", function(object) object@visible)
setMethod("color", "Interpolation", function(object) object@color)
setMethod("persistent", "Interpolation", function(object) object@persistent)


#--- replacement methods
setReplaceMethod("rates", "Interpolation", function(object, value) {
	object@rates<-value
	object
})

setReplaceMethod("visible", "Interpolation", function(object, value) {
	object@visible<-value
	object
})

setReplaceMethod("color", "Interpolation", function(object, value) {
	object@color<-value
	object
})

setReplaceMethod("persistent", "Interpolation", function(object, value) {
	object@persistent<-value
	object
})


#--- interpolate
setMethod("interpolate", c("Interpolation", "MareyMap"),	function(object, map) {
	object@physicalPositions <- physicalPositions(map)[which(markerValidity(map))]
	object@rates <- rep(0, length(physicalPositions(map)))
	object
})


#--- userParam
setMethod("userParam", "Interpolation", function(object) {	
		sav <- InterpolationParam()
		paramName(sav) <- "Saved"
		paramDesc(sav) <- "   Indicate if the interpolation is to be be kept   \n\ when the map is saved to text file"
		paramType(sav) <- "logical"
		paramDefault(sav) <- T
		paramFun(sav) <- "persistent"
    
    vis <- InterpolationParam()
		paramName(vis) <- "Displayed"
		paramDesc(vis) <- "   Whether the line is visible on the plot or not   "
		paramType(vis) <- "logical"
		paramDefault(vis) <- T
		paramFun(vis) <- "visible"
		
		col <- InterpolationParam()
		paramName(col) <- "Line color"
		paramDesc(col) <- "          Color of the line          "
		paramType(col) <- "color"
		paramDefault(col) <- "#000000"
		paramFun(col) <- "color"
	
    list(sav, vis, col)
	}
)


#--- query
setMethod("query", c("Interpolation", "numeric"), function(object, pos) {0})  ## replace 'integer' by 'numeric' to correct a bug


#--- argList
setMethod("argList", "Interpolation", function(object) {
  paste("visible = ", object@visible, ", persistent = ", object@persistent,	", color = \"", object@color, "\"",	sep = "")
})


#--- plot	
setMethod("plot", c("MareyMap", "Interpolation"), function(x, y, ...) {
	cl <- match.call()
	cl$col <- y@color
	cl$ylim <- cl$ylimr
	cl$x <- physicalPositions(x)[which(markerValidity(x))] / 1000000
	cl$y <- rates(y)[which(!is.na(rates(y)))]
	eval(cl)
})


#--- plotModel
setMethod("plotModel", "Interpolation", function(object, ...) {
	par(col = "white", col.lab = "white", col.axis = "white", tck = 0)
	plot(0, 0)
})


#--- plotRate
setMethod("plotRate", "Interpolation", function(object, ...) {
	cl <- match.call()
  plot(physicalPositions(object) / 1000000, rates(object)[which(!is.na(rates(object)))], col = object@color, ...)
})

