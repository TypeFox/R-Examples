`generate.ith.superclass` <- function(			#	specifying
	R,						#		a) the number of residual classes to bin into two superclasses
	superclass.index)				#		b) index of ith superclass composition of interest
{							#	the contents of the ith superclass is generated
	superclass.left <- numeric(R)
	superclass.left[1] <- 1						#first factor always present
	superclass.index <- superclass.index - 1			#use superclass.index as temporary variable
	for (superclass.left.index in 2:R) {
		superclass.left[superclass.left.index] <- superclass.index %% 2
		superclass.index <- floor(superclass.index / 2)
	}

	return(superclass.left)
}

