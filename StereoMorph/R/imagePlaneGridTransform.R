imagePlaneGridTransform <- function(p, nx, ny){

	# GET CORNERS FROM FIRST 8 INPUT PARAMETERS
	corners <- matrix(p[1:8], nrow=4, ncol=2, byrow=TRUE)

	# GET INTIAL AND FINAL ROW SLOPES
	rm1_m <- (corners[2, 2] - corners[1, 2])/(corners[2, 1] - corners[1, 1])
	rm2_m <- (corners[3, 2] - corners[4, 2])/(corners[3, 1] - corners[4, 1])

	# GET INTIAL AND FINAL COLUMN SLOPES
	cm1_m <- (corners[4, 2] - corners[1, 2])/(corners[4, 1] - corners[1, 1])
	cm2_m <- (corners[3, 2] - corners[2, 2])/(corners[3, 1] - corners[2, 1])
	
	# FIND POINTS BETWEEN CORNERS EITHER EVENLY SPACED OR AT LINEARLY CHANGING INTERVALS (CONSTANT 2ND DERIVATIVE)
	r_pos_x <- quadraticPointsOnInterval(t1=corners[1, 1], t2=corners[4, 1], n=ny, a=p[9])
	r_pos <- cbind(r_pos_x, r_pos_x*cm1_m + corners[4, 2] - cm1_m*corners[4, 1])
	c_pos_x <- quadraticPointsOnInterval(t1=corners[1, 1], t2=corners[2, 1], n=nx, a=p[10])
	c_pos <- cbind(c_pos_x, c_pos_x*rm1_m + corners[2, 2] - rm1_m*corners[2, 1])

	# SET SLOPE GRADIENTS FOR ROWS
	r_dygrad <- quadraticPointsOnInterval(t1=corners[2, 2] - corners[1, 2], t2=corners[3, 2] - corners[4, 2], n=ny, a=p[11])
	r_dxgrad <- quadraticPointsOnInterval(t1=corners[2, 1] - corners[1, 1], t2=corners[3, 1] - corners[4, 1], n=ny, a=p[11])

	# SET SLOPE GRADIENTS FOR COLS
	c_dygrad <- quadraticPointsOnInterval(t1=corners[4, 2] - corners[1, 2], t2=corners[3, 2] - corners[2, 2], n=nx, a=p[12])
	c_dxgrad <- quadraticPointsOnInterval(t1=corners[4, 1] - corners[1, 1], t2=corners[3, 1] - corners[2, 1], n=nx, a=p[12])

	# FILL GRID POINT PAIRWISE MATRICES WITH SLOPES
	rmm <- matrix(r_dygrad/r_dxgrad, nrow=nx, ncol=ny, byrow=T)
	cmm <- matrix(c_dygrad/c_dxgrad, nrow=nx, ncol=ny, byrow=F)

	# FILL GRID POINT PAIRWISE MATRICES WITH INTERCEPTS
	cbm <- matrix(c_pos[, 2] - cmm[, 1]*c_pos[, 1], nrow=nx, ncol=ny, byrow=F)
	rbm <- matrix(r_pos[, 2] - rmm[1, ]*r_pos[, 1], nrow=nx, ncol=ny, byrow=T)

	# WHAT IF M IS INF??
	#if(is.null(m1)) return(c(b1, m2*b1+b2))
	#if(is.null(m2)) return(c(b2, m1*b2+b1))

	# GET INTERSECTIONS BY PERFORMING CALCULATION ON PAIRWISE MATRICES
	x <- (cbm - rbm) / (rmm - cmm)
	y <- rmm * x + rbm

	# MAKE MATRIX
	grid <- cbind(c(x), c(y))
	
	# RETURN AS GRID
	grid
}