qsim.norm <-
function (p, mu, sig) {
	return(mu + qnorm(p) * sig * sqrt(mu^3 * (1-mu)^3))
}
