`simulateVAR` <-
function(R, T, P, v, perc) {

	## Create D matrix with perc% of possible edges non-null
	## Non-null edges are either U(0.2, 1) or U(-1, -0.2) with equal probability
	D <- matrix(0, nrow = P, ncol = P)
	index <- expand.grid(seq(1:P),seq(1:P))
	selected.index <- sample(seq(1:(P*P)), ceiling(0.10 * P * P))
	selected.edges <- index[selected.index,]
	for(edge in 1:ceiling(0.10 * P * P)) {
		tmp <- runif(1)
		if(tmp > 0.5) {
			D[selected.edges[edge,1], selected.edges[edge,2]] <-
			runif(1, min = 0.2, max = 1)
		}
		else {
			D[selected.edges[edge,1], selected.edges[edge,2]] <-
			runif(1, min = -1, max = -0.2)
		}
	}
	Dtrue <- abs(sign(D))

	## Simulate data
	y <- vector("list", R)
	for(r in 1:R) {
		y[[r]] <- matrix(NA, nrow = P, ncol = T)
		y[[r]][,1] <- rnorm(P, mean = 0, sd = sqrt(v^(-1)))
		for(t in 2:T) {
			y[[r]][,t] <- D %*% y[[r]][, t - 1] +
				as.matrix(rnorm(P,mean = 0, sd = sqrt(v^(-1))))
		}
	}
	return(list("Dtrue" = Dtrue, "y" = y))
}

