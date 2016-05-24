intersectCirclePlane <- function(circle, P, N){
	# http://mathforum.org/library/drmath/view/69136.html
	# http://mathforum.org/library/drmath/view/65138.html
	# P is a point in the plane
	# N is the normal vector of the plane

	# EQUATIONS CAN BE CONSOLIDATED TO JUST TWO SOLUTIONS
	# NOT NECESSARY TO SELECT CORRECT SOLUTIONS FROM FOUR
	# SEE angleOnCircleFromPoint

	# MAKE SURE N IS A UNIT VECTOR
	N <- uvector(N)

	# FIND COEFFICIENTS
	a <- sum(N*(circle$V*circle$R))
	b <- sum(N*(circle$U*circle$R))
	c <- -sum(N*(circle$C - P))

	alpha <- acos(a / sqrt(a^2 + b^2))
	t <- c(
		asin(c / sqrt(a^2 + b^2)) - alpha, 
		asin(c / sqrt(a^2 + b^2)) + alpha
		)

	a <- sum(-N*(circle$V*circle$R))
	b <- sum(-N*(circle$U*circle$R))
	c <- -sum(-N*(circle$C - P))

	alpha <- acos(a / sqrt(a^2 + b^2))
	t <- c(t,
		asin(c / sqrt(a^2 + b^2)) - alpha, 
		asin(c / sqrt(a^2 + b^2)) + alpha
		)

	# POINTS ON CIRCLE
	p <- circlePoint(circle, t)

	# FIND POINT IN PLANE
	t <- t[abs(distPointToPlane(p, N, P)) < 1e-9]

	if(length(t) == 1) t <- c(t,t)

	# FIND ORTHOGONAL PROJECTION OF CIRCLE CENTER IN PLANE
	t
}