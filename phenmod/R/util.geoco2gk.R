util.geoco2gk <- function(x,y, meridian=4) {
	sy <- meridian
	rho <- 180 / pi
	e2 <- 0.0067192188
	c <- 6398786.849

	for(i in 1:length(x)) {
    		bf <- as.numeric(y[i]) / rho
    		g <- 111120.61962 * y[i] -15988.63853 * 
			sin(2 * bf) + 16.72995 * 
			sin(4 * bf) - 0.02178 * 
			sin(6 * bf) +0.00003 * 
			sin(8 * bf)
    		co <- cos(bf)
    		g2 <- e2 * (co^2)
    		g1 <- c / sqrt(1+g2)
    		t_mod <- sin(bf) / cos(bf)
    		dl <- x[i] - sy * 3
    		fa <- co * dl / rho
    		y[i] <- g + (fa^2) * t_mod * 
			g1 / 2 + (fa^4) * t_mod * 
			g1 * (5 - (t_mod^2) + 9 * g2) / 24
    		rm_mod <- fa * g1 + (fa^3) * g1 * 
			(1 - (t_mod^2) + g2) / 6 + (fa^5) * 
			g1 * (5 - 18 * (t_mod^6)) / 120
    		x[i] <- rm_mod + sy * 1000000 + 500000
	}

 	result <- cbind(x,y)
 	return(result)
}

