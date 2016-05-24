FDR <-
function(test_statistics, t_0) {
	m2 <- mean(test_statistics^2)
	m4 <- mean(test_statistics^4)
	mu <- sqrt((3 - (6*m2) + m4)/(m2 - 1))
	theta <- ((m2 - 1)^2)/(3 - (6*m2) + m4)
	return((1 - theta) * pnorm(-t_0) / (((1 - theta) * pnorm(-t_0)) + (theta * (pnorm(mu - t_0) + pnorm(-mu - t_0)) / 2)))
	}
