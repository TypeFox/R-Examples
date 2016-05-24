# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MMLoess.
#
# MMLoess is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MMLoess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MMLoess; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


#--- class definition
setClass("MMLoess",
  contains = "Interpolation",
  representation(
    name = "character",
    span = "numeric",
  	degree = "integer",
    model = "ANY"
  ),
  prototype(
    name = "loess",
    span = 0.1,
		degree = as.integer(2),
    model = NULL
  )
)


#--- constructors
setGeneric("MMLoess", function(val) standardGeneric("MMLoess"))
setMethod("MMLoess", signature(val = "missing"), function(val) {
  object <- new("MMLoess")
  object@visible <- T
  object@color <- "#000000"
  object
})


#--- accessors
setMethod("name", "Interpolation", function(object) object@name)

setGeneric("span", function(object) standardGeneric("span"))
setMethod("span", "MMLoess", function(object) object@span)

setGeneric("degree", function(object) standardGeneric("degree"))
setMethod("degree", "MMLoess", function(object) object@degree)


#--- replacement methods
setReplaceMethod("name", "Interpolation",	function(object, value) {
	    object@name<-value
	    object
  })

setGeneric("span<-", function(object, value) standardGeneric("span<-"))
setReplaceMethod("span", "MMLoess", function(object, value) {
  if(value > 0)
	  object@span <- value
  object
})


setGeneric("degree<-", function(object, value) standardGeneric("degree<-"))
setReplaceMethod("degree","MMLoess", function(object, value) {
  if(value > 0 && value < 3)
  	object@degree <- as.integer(value)
  object
})


#--- query
setMethod("query", c("MMLoess", "numeric"), function(object, pos) {
  sapply(pos, function(pp) {
    posplus <- pp + 1
    posmoins <- pp - 1
    pred <- predict(object@model, newdata = c(posmoins, posplus))        
    round((pred[[2]] - pred[[1]]) / ((posplus - posmoins) / 1000000), digits = 2)
	})
})


#--- interpolate
setMethod("interpolate", signature(object = "MMLoess", map = "MareyMap"), function(object, map) {
  valgen <- geneticDistances(map)[which(markerValidity(map))]
  valphys <- physicalPositions(map)[which(markerValidity(map))]
  object@physicalPositions <- valphys
  object@model <- loess(valgen ~ valphys, span = object@span, degree = object@degree, na.action = na.exclude)
  ## physical positions + 1
  pp1 <- physicalPositions(map) + 1
  ## physical positions - 1
  pm1 <- physicalPositions(map)- 1
  ## predicted genetical positions @ phys + 1
  gp1 <- predict(object@model, newdata = pp1)
  ## predicted genetical positions @ phys - 1
  gm1 <- predict(object@model, newdata = pm1)
  object@rates <- mapply(function(xa, xb, ya, yb) {round((yb - ya) / ((xb - xa) / 1000000), 2)}, pm1, pp1, gm1, gp1)
  object
})


#--- plot
setMethod("plot", c("MareyMap", "MMLoess"), function(x, y, ...) {
  plotModel(y, ...)
  par(new = T)
  plotRates(y, ...)
})


#--- plotModel
setMethod("plotModel", "MMLoess", function(object, ...) {
  expli <- seq(from = min(object@physicalPositions), to = max(object@physicalPositions), length.out = 100)
  plot(expli / 1000000, predict(object@model, newdata = expli), col = object@color, type = "l", ...)
})

#--- plotRate
setMethod("plotRate", "MMLoess", function(object, ...) {
  cl <- match.call()
  cl[[1]] <- as.name("plot")
  cl$object <- NULL
  cl$col <- object@color
  cl$type <- "l"
  expli <- seq(from = min(object@physicalPositions), to = max(object@physicalPositions), length.out = 100)
  cl$x <- expli / 1000000
  ## physical positions + 1
  pp1 <- expli + 1
  ## physical positions - 1
  pm1 <- expli - 1
  ## predicted genetical positions @ phys + 1
  gp1 <- predict(object@model, newdata = pp1)
  ## predicted genetical positions @ phys - 1
  gm1 <- predict(object@model, newdata = pm1)
  cl$y <- mapply(function(xa, xb, ya, yb){(yb - ya) / ((xb - xa) / 1000000)}, pm1, pp1, gm1, gp1)
  eval(cl)
})


#--- createOrder
setMethod("createOrder", "MMLoess", function(object) {
  paste("new(\"MMLoess\", name = \"", object@name, "\", ", argList(object), ", span = ", object@span, ", degree = as.integer(", object@degree, "))", sep = "")
})


#--- userParam
setMethod("userParam", "MMLoess", function(object) {
  nam <- InterpolationParam()
  paramName(nam) <- "Name"
	paramDesc(nam) <- "   The name of the interpolation   "
	paramType(nam) <- "character"
	paramDefault(nam) <- "loess"
	paramMin(nam) <- 1
	paramMax(nam) <- NULL	
	paramFun(nam) <- "name"

	span <- InterpolationParam()
  paramName(span) <- "Span"
  paramDesc(span) <- "   Span controls the number of surrounding markers   \n\ that are to be taken into account while fitting\n\ a local polynomial. It is defined as a percentage\n\ of the complete set of markers. The higher it is,\n\ the smoother the curve gets.\n\ \n\ See \"Using MareyMap\" for more informations."
  paramType(span) <- "numeric"
  paramDefault(span) <- 0.1
  paramMin(span) <- 0
  paramMax(span) <- 1    
  paramFun(span) <- "span"
        
  deg <- InterpolationParam()
  paramName(deg) <- "Degree"
  paramDesc(deg) <- "   Degree of the local fitted polynomials, up to 2   "
  paramType(deg) <- "integer"
  paramDefault(deg) <- 2
  paramMin(deg) <- 0
  paramMax(deg) <- 2    
  paramFun(deg) <- "degree"
  
  c(list(nam), userParam(as(object, "Interpolation")), list(span, deg))
})

registerInterpolationMethod("Loess based method", "MMLoess")
