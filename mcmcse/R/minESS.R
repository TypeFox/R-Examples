minESS <- function(p, alpha = .05, eps = .05)
{

	crit <- qchisq(1-alpha, p)
	foo <- 2/p
	logminESS <- foo*log(2) + log(pi) - foo*log(p) - foo*lgamma(p/2) - 2*log(eps) + log(crit)
	return(round(exp(logminESS)))

}