###
### CROSS.R  Vector product
###

cross <- function(x, y) {
	if (!is.numeric(x) || !is.numeric(y))
		stop("Arguments 'x' and 'y' must be numeric vectors or matrices.")

	if (is.vector(x) && is.vector(y)) {
		if (length(x) == length(y) && length(x) == 3) {
			xxy <- c(x[2]*y[3] - x[3]*y[2],
			         x[3]*y[1] - x[1]*y[3],
    		         x[1]*y[2] - x[2]*y[1])
		} else {
			stop("Vectors 'x' and 'y' must be both of length 3.")
		}
	} else {
		if (is.matrix(x) && is.matrix(y)) {
			if (all(dim(x) == dim(y))) {
				if (ncol(x) == 3) {
                    xxy <- cbind(x[, 2]*y[, 3] - x[, 3]*y[, 2],
        			             x[, 3]*y[, 1] - x[, 1]*y[, 3],
            		             x[,1 ]*y[, 2] - x[, 2]*y[, 1])
				} else {
					if (nrow(x) == 3) {
                        xxy <- rbind(x[2, ]*y[3, ] - x[3, ]*y[2, ],
            			             x[3, ]*y[1, ] - x[1, ]*y[3, ],
                		             x[1, ]*y[2, ] - x[2, ]*y[1, ])
					} else {
						stop("'x', 'y' must have one dimension of length 3.")
					}
				}
			} else {
				stop("Matrices 'x' and 'y' must be of same size.")
			}
		} else {
			if (is.vector(x) && is.matrix(y) ||
				is.matrix(x) && is.vector(y)) {
			stop("Arguments 'x', 'y' must be vectors/matrices of same size.")
				}
		}
	}
	return(xxy)
}
