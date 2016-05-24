soilwater3 <- 
function(x, theta_R, a1, p1, a2, p2)
	theta_R + a1 * exp(-x/p1) + a2 * exp(-x/p2)
