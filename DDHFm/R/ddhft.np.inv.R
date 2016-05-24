"ddhft.np.inv" <-
function(hft.obj) {

#	initialisation

	hft <- hft.obj$hft
	factors <- hft.obj$factors
	n <- length(hft)
	nhalf <- n/2
	J <- logb(n, 2)
	sm <- rep(0, nhalf)
	det <- sm
	
#	decomposition part

	for (i in 1:J) {
		sm[1:nhalf] <- (hft[2 * (1:nhalf) - 1] + hft[2 * (1:nhalf)])/2
		det[1:nhalf] <- (hft[2 * (1:nhalf) - 1] - hft[2 * (1:nhalf)])/2
		hft[1:nhalf] <- sm[1:nhalf]
		hft[(nhalf+1):n] <- det[1:nhalf]
		n <- n/2
		nhalf <- n/2
	}

	nhalf <- 1
	n <- 2

#	reconstruction part

	for (i in 1:J) {
		sm[1:nhalf] <- hft[1:nhalf]
		det[1:nhalf] <- hft[(nhalf+1):n]
		v <- factors[(nhalf+1):n]
		det[1:nhalf] <- det[1:nhalf] * v
		hft[2 * (1:nhalf) - 1] <- sm[1:nhalf] + det[1:nhalf]
		hft[2 * (1:nhalf)] <- sm[1:nhalf] - det[1:nhalf]
#		hft[1:n][hft[1:n] < 0] <- 0
		nhalf <- n
		n <- 2 * n
	}
	return(hft)
}

