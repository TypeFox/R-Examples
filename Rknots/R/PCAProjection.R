##Script generated in:
# 2011
# 5:12:17 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

PCAProjection <- function (points3D) {
	nrow <- nrow(points3D)
	ncol <- ncol(points3D)
	centroid <- apply(points3D, 2, mean)
	for (i in 1 : nrow) 
		points3D[i, ] <- points3D[i, ] - centroid
	covM = t(points3D) %*% points3D
	eigenvectors <- t(eigen(covM)$vectors)
	if (det(eigenvectors) < 0) 
		eigenvectors[3, ] = -1 * eigenvectors[3, ]
	points3D.rot <- points3D %*% t(eigenvectors)
	return(points3D.rot)
}



