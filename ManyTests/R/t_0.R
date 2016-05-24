t_0 <-
function(test_statistics) {
	m2 <- mean(test_statistics^2)
	m4 <- mean(test_statistics^4)
	mu <- sqrt((3 - (6*m2) + m4)/(m2 - 1))
	theta <- ((m2 - 1)^2)/(3 - (6*m2) + m4)
	return((1/mu) * log((((1 - theta)/theta) * exp((mu^2)/2)) + sqrt(((((1 - theta)/theta)^2) * exp(mu^2)) - 1)))
	}
