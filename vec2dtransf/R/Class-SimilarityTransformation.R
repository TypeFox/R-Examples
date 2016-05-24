### SimilarityTransformation class ###
## Similarity transformations can be written as:
##      x’ = ax + by + c
##      y’ = ay - bx + d
setClass("SimilarityTransformation",
	representation(),
	contains = "Cartesian2DCoordinateTransformation",
	validity = function(object) {
		if (length(object@parameters) == 0){
			if (nrow(object@controlPoints) < 2)
				stop("At least 2 control points (rows in the data.frame 'controlPoints') are required for the similarity transformation")
		}
		else if (length(object@parameters) != 4){ 
			stop("Similarity transformations require 4 parameters!")
		}
		return(TRUE)
	}
)
