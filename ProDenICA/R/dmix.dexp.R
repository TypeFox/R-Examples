dmix.dexp <-
function(x, means = 0)
{
	f1 <- x * 0
	nm <- length(means)
	for(m in means) {
		xx <- x + m
		f1 <- f1 + (exp( - xx) * (xx > 0) + exp(xx) * (xx < 0))/
			nm
	}
	f1
}

