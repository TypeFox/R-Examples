# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MareyMap.
#
# MMSpline3 is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MMSpline3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MMSpline3; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


#--- class definition
setClass("MMSpline3",	
  contains = "Interpolation",
  representation(
    name = "character",
  	type = "character",
    gcv = "logical",
    df = "numeric",
    spar = "numeric",
    model = "ANY"
  ),
	prototype(
    name = "spline",
    type = "cross-validation",
    df = 10,
    spar = 0.1,
    gcv = FALSE,
    model = NULL
  )
)


##--- constructors
setGeneric("MMSpline3",	function(val) standardGeneric("MMSpline3"))
setMethod("MMSpline3", signature(val = "missing"), function(val) {
	object <- new("MMSpline3")
	object@visible <- T
	object@color <- "#000000"
	object
})


##--- accessors
setMethod("name", "Interpolation", function(object) object@name)

setGeneric("type", function(object) standardGeneric("type"))
setMethod("type", "MMSpline3", function(object) object@type)

setGeneric("spar", function(object) standardGeneric("spar"))
setMethod("spar","MMSpline3",	function(object) object@spar)

setGeneric("df", function(object) standardGeneric("df"))
setMethod("df","MMSpline3",	function(object) object@df)

setGeneric("gcv", function(object) standardGeneric("gcv"))
setMethod("gcv", "MMSpline3", function(object) object@gcv)


##--- replacement methods
setReplaceMethod("name", "Interpolation",	function(object, value) {
	    object@name<-value
	    object
  })

setGeneric("type<-", function(object, value) standardGeneric("type<-"))
setReplaceMethod("type", c("MMSpline3", "character"), function(object, value) {
	if(value == "spar" || value == "degree of freedom" || value == "cross-validation")
  	object@type <- value
  object
})


setGeneric("spar<-", function(object, value) standardGeneric("spar<-"))
setReplaceMethod("spar", c("MMSpline3", "numeric"), function(object, value) {
	if(value > 0)
		object@spar <- value
	object
})


setGeneric("df<-", function(object, value) standardGeneric("df<-"))
setReplaceMethod("df", c("MMSpline3", "numeric"),	function(object, value) {
	if(value > 0)
		object@df <- value
		object
})


setGeneric("gcv<-",	function(object, value) standardGeneric("gcv<-"))
setReplaceMethod("gcv", c("MMSpline3", "logical"), function(object, value) {
	object@gcv <- value
	object
})


##--- query
setMethod("query", c("MMSpline3", "numeric"), function(object, pos) {
  sapply(pos, function(pp) {
    posplus <- pp + 1
    posmoins <- pp - 1
    pred <- predict(object@model, c(posmoins, posplus))$y
    round((pred[[2]] - pred[[1]]) / ((posplus - posmoins) / 1000000), digits = 2)
  })
})


##--- interpolate
setMethod("interpolate", signature(object = "MMSpline3", map = "MareyMap"),	function(object, map) {
	valgen <- geneticDistances(map)[which(markerValidity(map))]
	valphys <- physicalPositions(map)[which(markerValidity(map))]
	object@physicalPositions <- valphys
  cl <- call("smooth.spline", x = as.vector(valphys), y = as.vector(valgen))
  if(object@type == "cross-validation") {
  	if(object@gcv)
    	cl$cv <- FALSE
    else
      cl$cv <- TRUE
  } else if(object@type == "spar") {
    cl$spar <- object@spar
  } else if(object@type == "degree of freedom") {
    cl$df <- object@df
  }
	object@model <- eval(cl)
  object@rates <- query(object, physicalPositions(map))
  object
})


##--- plot
setMethod("plot", c("MareyMap", "MMSpline3"), function(x, y, ...) {
	plotModel(y, ...)
	par(new = T)
	plotRates(y, ...)
})


##--- plotModel  # Method called for plotting the interpolation model (if any)
setMethod("plotModel", "MMSpline3", function(object, ...) {
	cl <- match.call()
	cl[[1]] <- as.name("plot")
  cl$object <- NULL
	cl$col <- object@color
	cl$type <- "l"
	expli <- seq(from = min(object@physicalPositions), to = max(object@physicalPositions), length.out = 100)
	cl$x <- expli / 1000000
	cl$y <- predict(object@model, expli)$y
	eval(cl)
})


##--- plotRate # method called to plot the recombination rates
setMethod("plotRate", "MMSpline3", function(object, ...) {
	cl <- match.call()
	cl[[1]] <- as.name("plot")
	cl$object <- NULL
	cl$col <- object@color
	cl$type <- "l"
	expli <- seq(from = min(object@physicalPositions), to = max(object@physicalPositions), length.out = 100)
	cl$x <- expli / 1000000
  cl$y <- query(object, expli)
  invisible(cl)
  eval(cl)
})


##--- createOrder
setMethod("createOrder", "MMSpline3",	function(object) {
  paste("new(\"MMSpline3\", name = \"", object@name, "\", ", argList(object), ", spar = ", object@spar, ", df = ", object@df, ")", sep = "")
})


##--- userParam
setMethod("userParam", "MMSpline3",	function(object) {
  nam <- InterpolationParam()
  paramName(nam) <- "Name"
  paramDesc(nam) <- "   The name of the interpolation   "
  paramType(nam) <- "character"
  paramDefault(nam) <- "spline"
  paramMin(nam) <- 1
  paramMax(nam) <- NULL	
  paramFun(nam) <- "name"
    
	type <- InterpolationParam()
	paramName(type) <- "Type"
	paramDesc(type) <- "   There are three different ways of estimating the   \n\ parameters of locally fitted polynomials. 
If you choose 'spar' here, you will have to set\n\ the value for the parameter 'spar' below.
If you select 'degree of freedom', you will have to set\n\ the value of the parameter 'degree of freedom'. 
If you choose 'cross-validation', you have to\n\ choose whether or not to use 'generalised \n\ cross-validation'.\n\
Using 'cross-validation' with generalised cross-\n\ validation set off is generally a good choice on\n\ a clean map."
	paramType(type) <- "list"
  paramValues(type) <- c("spar", "degree of freedom", "cross-validation")
  paramDefault(type) <- "cross-validation"
	paramFun(type) <- "type"

	spar <- InterpolationParam()
	paramName(spar) <- "Spar"
	paramDesc(spar) <- "   Spar is the smoothing parameter, similar to span in Loess.   
The higher it is, the smoother the curve gets."
	paramType(spar) <- "numeric"
	paramDefault(spar) <- 0.1
	paramMin(spar) <- 0.0000001
	paramMax(spar) <- 1
	paramFun(spar) <- "spar"
	
	df <- InterpolationParam()
	paramName(df) <- "Degree of freedom"
	paramDesc(df) <- "   Controls the amount of smoothing by setting the   \n\ degree of freedom, which corresponds to the trace\n\ of the smoothing matrix. It is not intuitively\n\ obvious how to choose a value for this parameter\n\ and is often more convinient to rely on spar or\n\ cross-validation."
	paramType(df) <- "integer"
	paramDefault(df) <- 2
	paramMin(df) <- 0
	paramMax(df) <- 100	
	paramFun(df) <- "df"

  gcv <- InterpolationParam()
	paramName(gcv) <- "Generalised cross-validation"
	paramDesc(gcv) <- "   Whether or not Generalised Cross Validation shall   \n\ be used instead of usual cross validation. The\n\ 'generalized' cross-validation method should be\n\ used when there are duplicated points in 'x'."
	paramType(gcv) <- "logical"
	paramDefault(gcv) <- FALSE
	paramFun(gcv) <- "gcv"
  
	c(list(nam), userParam(as(object, "Interpolation")), list(type, spar, df, gcv))
})

registerInterpolationMethod("Cubic splines", "MMSpline3")
