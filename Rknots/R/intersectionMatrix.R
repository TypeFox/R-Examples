# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

intersectionMatrix <- function (points3D, ends = c()) 
{
	if (missing(points3D)) 
		stop("intersectionMatrix: Argument 'points3D' missing, with no default\n")
	nedge <- nrow(points3D) - 1
	M <- matrix(0, nedge, nedge)
	int.list <- list()
	for (i in 1 : nedge) {
		int.list[[i]] <- edgeIntersections(points3D[1 : (i + 1), ], 1 : (i + 1), i)
	}
	n.int <- unlist(lapply(int.list, length)) / 2
	intersecting <- which(n.int != 0)
	hits <- length(intersecting)
	n.int <- n.int[intersecting]
	for (i in 1 : hits) {
		for (j in 1 : n.int[i]) {
			M[intersecting[i], int.list[[intersecting[i]]][j, 1]] <- 
					int.list[[intersecting[i]]][j, 2]
		}
	}
	M <- M - t(M)
	M[ends, ] <- M[ends, ] * 0
	M[, ends] <- M[, ends] * 0
	return(M)
}

