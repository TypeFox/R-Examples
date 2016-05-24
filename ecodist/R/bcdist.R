bcdist <- function(x, rmzero = FALSE)
{

# Calculates Bray-Curtis distance (1 - percent similarity) for the
# rows of a matrix.
# 
# If rmzero = T, all empty rows of the matrix will be removed before calculating
# distances. Otherwise the distance between empty rows will be set to 0.
#
# Sarah C. Goslee, November 1997
#
#
# To use this function, two files are necessary.
# This file is bcdist.q
# The other is bcdist.c (incorporated in ecodist.c).

	x <- as.matrix(x)
	if(rmzero == TRUE) {
		xsum <- apply(x, 1, sum)
		x <- x[xsum > 0,  ]
	}
	dist.v <- rep(0, (nrow(x) * (nrow(x) - 1))/2)	

   cresults <- .C("bcdist",
		as.double(as.vector(t(x))),
		as.integer(nrow(x)),
		as.integer(ncol(x)),
		dist.v = as.double(dist.v),
		PACKAGE = "ecodist")
	dist.v <- cresults$dist.v

	## give the results attributes similar to dist()
	attr(dist.v, "Size") <- nrow(x)
 	 attr(dist.v, "Labels") <- rownames(x)
	 attr(dist.v, "Diag") <- FALSE
	 attr(dist.v, "Upper") <- FALSE
    attr(dist.v, "method") <- "bray-curtis"
    class(dist.v) <- "dist"

	dist.v
}
