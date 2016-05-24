pointPlaneProj <- function(q, p, n){

	if(q[1] == "non-coincident"){
		warning("Input point non-coincident")
		return(c(NA, NA, NA))
	}

	# PROJECTION OF POINT Q ONTO PLANE DEFINED BY POINT P AND NORMAL VECTOR N
	q - sum((q - p) * n) * n
}