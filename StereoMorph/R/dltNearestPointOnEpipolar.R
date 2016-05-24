dltNearestPointOnEpipolar <- function(p1, p2, cal.coeff){

	# GET EPIPOLAR LINE
	epipolar_line <- dltEpipolarLine(p1, cal.coeff)

	if(is.vector(p2)){

		# GET DISTANCE FROM SECOND POINT TO EPIPOLAR LINE
		edist <- dltEpipolarDistance(p1=p1, p2=p2, cal.coeff=cal.coeff)

		# GET CLOSEST POINT ON EPIPOLAR LINE
		matching.pt <- orthogonalProjectionToLine(p=p2, l1=epipolar_line$l1, l2=epipolar_line$l2)

		return(list(matching.pt=matching.pt, min.idx=1, p2.dist=edist))
	}

	# FIND DISTANCES FROM EPIPOLAR LINE TO SECOND SET OF POINTS
	edist <- dltEpipolarDistance(p1=p1, p2=p2, cal.coeff=cal.coeff)

	# GET SECOND POINT AT MINIMUM DISTANCE
	min.idx <- which.min(edist)

	# GET CLOSEST POINT ON EPIPOLAR LINE
	matching.pt <- orthogonalProjectionToLine(p=p2[min.idx, ], l1=epipolar_line$l1, l2=epipolar_line$l2)

	# GET MINIMUM DISTANCE FROM SECOND POINT SET
	p2.dist <- edist[min.idx]

	# DISTANCE FROM EPIPOLAR LINE IS ZERO
	epipolar_dist <- 0

	list(matching.pt=matching.pt, min.idx=min.idx, p2.dist=p2.dist)
}