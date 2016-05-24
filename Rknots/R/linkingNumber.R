# Script comments and history
# 2011
# Feb 11, 2011
# 8:39:54 AM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

linkingNumber <- function(points3D, ends, M = c())
{
	nedge <- nrow(points3D) - 1
	if (is.null(M)) {
		M <- intersectionMatrix(points3D)		
	}
	subM <- M[1 : (ends[1] - 1), (ends[1] + 1) : nedge]
	ones <- which(subM != 0, arr.ind = TRUE)
	ones[, 2] <- ones[, 2] + ends[1]
	skein.signs <- apply(ones, 1, skeinSign, points3D = points3D)
	lk <- sum(skein.signs) / 2
	return (lk)	
}
