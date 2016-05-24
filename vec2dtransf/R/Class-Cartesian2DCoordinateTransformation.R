setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass("Cartesian2DCoordinateTransformation",
	representation(controlPoints = "data.frame", parameters="numeric", 
        residuals="matrix", rmse="numericOrNULL", "VIRTUAL"),
	prototype(controlPoints=data.frame(),parameters=numeric(),
        residuals=matrix(nrow=0,ncol=0),rmse=NULL),
	validity = function(object) {
		if (!inherits(object@controlPoints, "data.frame"))
			stop("controlPoints should be of class data.frame")
		if (ncol(object@controlPoints) == 0 && length(object@parameters) == 0 )
			stop("Either 'controlPoints' or 'parameters' must be provided!")
		if ( length(object@parameters) == 0 ){
			if (ncol(object@controlPoints) != 4)
				stop("'controlPoints' must have four (4) columns: X source, Y source, X target, Y target")
		}
		return(TRUE)
	}
)
