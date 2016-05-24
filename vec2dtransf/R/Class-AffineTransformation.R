### AffineTransformation class ###
## Affine transformations can be written as:
##      x’ = ax + by + c
##      y’ = dx + ey + f
setClass("AffineTransformation",
	representation(),
	contains = "Cartesian2DCoordinateTransformation",
	validity = function(object) {
		if (length(object@parameters) == 0){
			if (nrow(object@controlPoints) < 3)
				stop("At least 3 control points (rows in the data.frame 'controlPoints') are required for the affine transformation")
		}
		else if (length(object@parameters) != 6){ 
			stop("Affine transformations require 6 parameters!")
		}
		return(TRUE)
	}
)

