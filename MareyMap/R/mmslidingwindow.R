# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MMSlidingWindow.
#
# MMSlidingWindow is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MMSlidingWindow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MMSlidingWindow; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


#--- class definition
setClass("MMSlidingWindow",
  contains = "Interpolation",
  representation(
    name = "character",
    size = "integer",
    shift = "integer",
    threshold = "integer",
    results = "vector",
    geneticalDistances = "vector"
  ),
  prototype(
    name = "slidingwindow",
    size = as.integer(1000000),
    shift = as.integer(500000),
    threshold = as.integer(8),
    results = vector(),
    geneticalDistances = vector()
  )
)


#--- constructors
setGeneric("MMSlidingWindow", function(val) standardGeneric("MMSlidingWindow"))
setMethod("MMSlidingWindow", signature(val = "missing"), function(val) {
  object <- new("MMSlidingWindow")
  object@visible <- T
  object@color <- "#000000"
  object
})


#--- accessors
setMethod("name", "Interpolation", function(object) object@name)

setGeneric("size", function(object) standardGeneric("size"))
setMethod("size", "MMSlidingWindow", function(object) object@size)

setGeneric("shift", function(object) standardGeneric("shift"))
setMethod("shift", "MMSlidingWindow", function(object) object@shift)

setGeneric("threshold", function(object) standardGeneric("threshold"))
setMethod("threshold", "MMSlidingWindow", function(object) object@threshold)


#--- replacement methods
setReplaceMethod("name", "Interpolation",	function(object, value) {
	    object@name<-value
	    object
  })

setGeneric("size<-", function(object, value) standardGeneric("size<-"))
setReplaceMethod("size", "MMSlidingWindow", function(object, value) {
  if(value > 0)
	  object@size <- as.integer(value)
  object
})

setGeneric("shift<-", function(object, value) standardGeneric("shift<-"))
setReplaceMethod("shift","MMSlidingWindow", function(object, value) {
  if(value > 0)
    object@shift <- as.integer(value)
  object
})

setGeneric("threshold<-", function(object, value) standardGeneric("threshold<-"))
setReplaceMethod("threshold", "MMSlidingWindow", function(object, value) {
  if(value > 0)
    object@threshold <- as.integer(value)
  object
})


#--- query
setMethod("query", c("MMSlidingWindow", "numeric"), function(object, pos) {
  sapply(pos, function(pp) {
    if(pp < min(object@physicalPositions, na.rm = TRUE) || pp > max(object@physicalPositions, na.rm = TRUE))
      return(NA) 
    btm <- pp - object@size / 2
    top <- pp + object@size / 2
    idx <- which((object@physicalPositions >= btm) & (object@physicalPositions <= top))
    if(length(idx) < object@threshold)
      return(NA)
    mkrphys <- object@physicalPositions[idx]
    mkrgen <- object@geneticalDistances[idx]
    lm1 <- lm(mkrgen ~ mkrphys)
    round(lm1$coef[[2]] * 1000000, digits = 2)
  })
})


#--- queryModel -- non exported
setGeneric("queryModel", function(object, pos) standardGeneric("queryModel"))
setMethod("queryModel", c("MMSlidingWindow","numeric"), function(object, pos) {
  sapply(pos, function(pp) {
    if(pp < min(object@physicalPositions, na.rm = TRUE) || pp > max(object@physicalPositions, na.rm = TRUE))
      return(NA) 
    btm <- pp - object@size / 2
    top <- pp + object@size / 2
    idx <- which((object@physicalPositions >= btm) & (object@physicalPositions <= top))
    if(length(idx) < object@threshold)
      return(NA)
    mkrphys <- object@physicalPositions[idx]
    mkrgen <- object@geneticalDistances[idx]
    lm1 <- lm(mkrgen ~ mkrphys)
    lm1$coef[[1]] + lm1$coef[[2]] * pp
  })
})


#--- interpolate
setMethod("interpolate", signature(object = "MMSlidingWindow", map = "MareyMap"), function(object, map) {
  end <- max(physicalPositions(map))
  valgen <- geneticDistances(map)[which(markerValidity(map))]
  valphys <- physicalPositions(map)[which(markerValidity(map))]
  object@physicalPositions <- valphys
  object@geneticalDistances <- valgen
  object@rates <- query(object, physicalPositions(map))
  object
})


#--- plot
setMethod("plot", c("MareyMap", "MMSlidingWindow"), function(x, y, ...) {
	plotModel(y, ...)
  par(new = T)
  plotRates(y, ...)
})


#--- plotModel
setMethod("plotModel", "MMSlidingWindow", function(object, ...) {
  cl <- match.call()
  cl[[1]] <- as.name("plot")
  cl$object <- NULL
  cl$col <- object@color
  cl$type <- "l"
  expli <- seq(object@size / 2, max(object@physicalPositions), object@shift)
  cl$x <- expli / 1000000
  cl$y <- queryModel(object, expli)
  eval(cl)
})


#--- plotRate
setMethod("plotRate", "MMSlidingWindow", function(object, ...) {
  cl <- match.call()
  cl[[1]] <- as.name("plot")
  cl$object <- NULL
  cl$col <- object@color
  cl$type <- "l"
  expli <- seq(object@size / 2, max(object@physicalPositions), object@shift)
  cl$x <- expli / 1000000
  cl$y <- query(object, expli)
  eval(cl)
})


#--- createOrder
setMethod("createOrder", "MMSlidingWindow", function(object) {
	paste("new(\"MMSlidingWindow\", name = \"", object@name, "\", ", argList(object), ", size = as.integer(", object@size, "), shift = as.integer(", object@shift, "), threshold = as.integer(", object@threshold, "))", sep = "")
})


#--- userParam
setMethod("userParam", "MMSlidingWindow", function(object) {	
	nam <- InterpolationParam()
  paramName(nam) <- "Name"
  paramDesc(nam) <- "   The name of the interpolation   "
  paramType(nam) <- "character"
  paramDefault(nam) <- "slidingwindow"
  paramMin(nam) <- 1
  paramMax(nam) <- NULL	
  paramFun(nam) <- "name"

  size <- InterpolationParam()
  paramName(size) <- "Size"
  paramDesc(size) <- "   Size of the sliding window in base pair   "
  paramType(size) <- "integer"
  paramDefault(size) <- 3000000
  paramMin(size) <- 2
  paramFun(size) <- "size"

  shift <- InterpolationParam()
  paramName(shift) <- "Shift"
  paramDesc(shift) <- "   Shift between two consecutive windows in base pair   "
  paramType(shift) <- "integer"
  paramDefault(shift) <- 500000
  paramMin(shift) <- 1
  paramFun(shift) <- "shift"

  thr <- InterpolationParam()
  paramName(thr) <- "Threshold"
  paramDesc(thr) <- "   Minimum number of markers by window   "
  paramType(thr) <- "integer"
  paramDefault(thr) <- 5
  paramMin(thr) <- 2	
  paramFun(thr) <- "threshold"

  c(list(nam), userParam(as(object, "Interpolation")), list(size, shift, thr))
})

registerInterpolationMethod("Fixed size sliding window method", "MMSlidingWindow")
