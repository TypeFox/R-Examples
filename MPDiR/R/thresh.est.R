thresh.est <- function(p, obj) {
	cc <- coef(obj)
	m <- -cc[1]/cc[2]
	std <- 1/cc[2]
	qnorm(p, m, std)
	}