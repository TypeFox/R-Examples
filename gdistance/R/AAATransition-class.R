# Author: Jacob van Etten, jacobvanetten@yahoo.com
# International Rice Research Institute, IE University
# Date :  December 2010
# Version 1.0
# Licence GPL v3

setClass(Class="TransitionData",
		representation = representation(
			transitionMatrix = "sparseMatrix",
			transitionCells = "numeric",
			matrixValues = "character"
		),
    prototype(
      transitionMatrix = Matrix(0,1,1),
      transitionCells = 1,
      matrixValues = "conductance"
    ),
		validity = function(object){
			cond1 <- (nrow(transitionMatrix(object)) == ncol(transitionMatrix(object))) 
			cond2 <- (object@matrixValues == "resistance" | object@matrixValues == "conductance")
			cond3 <- length(transitionCells(object)) == object@transitionMatrix@Dim[1]
			cond <- cond1 & cond2 & cond3 
			return(cond)
	}
)

setClass(Class="TransitionLayer",
		contains = c("BasicRaster", "TransitionData"),
    prototype = prototype(
      rotated = FALSE,
      ncols= as.integer(1),
      nrows= as.integer(1),
      layernames=c(""),
      unit=c(""),
      z = list(),
      crs = CRS(as.character(NA)),
      transitionMatrix = Matrix(0,1,1),
      transitionCells = 1,
      matrixValues = "conductance"
      ),
		validity = function(object){
			cond1 <- (nrow(object@transitionMatrix) == ncol(object@transitionMatrix)) 
			cond2 <- (object@matrixValues == "resistance" | object@matrixValues == "conductance")
			cond3 <- length(transitionCells(object)) == object@transitionMatrix@Dim[1]
			cond <- cond1 & cond2 & cond3 
			return(cond)
	}
)

setClass ("TransitionStack",
	contains = "BasicRaster",
	representation (
			transition = "list"
	),
  prototype(
    rotated = FALSE,
    ncols= as.integer(1),
    nrows= as.integer(1),
    layernames=c(""),
    unit=c(""),
    z = list(),
    crs = CRS(as.character(NA)),
    transition=list(new("TransitionData"))
  ),
	validity = function(object) {return(TRUE)}
)

setClassUnion("Transition", c("TransitionLayer", "TransitionStack"))
