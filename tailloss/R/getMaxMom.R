getMaxMom <- function (x, theta = 0, cap = Inf, s, t = 1) {

  stopifnot(x >= 0,
  			is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0,
            is.scalar(cap), cap > 0)
            
	mom <- getMoments(x$Loss, theta = theta, cap = cap, maxmom = 2)
	lambda <- sum(x$Rate)
	EY <- (crossprod(x$Rate, mom)) / lambda
	ES <- lambda * t * EY

	# Approximate S with Gamma with same first two moments
	
	alphaS <- ES[1]^2 / ES[2]
	betaS <- ES[1] / ES[2]

	# find an interval that contains the optimal k for the moment bound of the Gamma approximation of S
	
	f <- function(k) lgamma(alphaS + k) - k * log(betaS) - k * log(s)
	
	oldk <-1
	k <- 2
	
	repeat {
		newk <- choose(k+1, 2) 	# newk grows following the triangular number sequence
		
		if (f(newk) > f(oldk))  break 	# since f is a convex function we stop when f starts to increase
		
		else oldk <- newk; k <- k + 1
		}
		
	maxmom <- newk

	maxmom
}