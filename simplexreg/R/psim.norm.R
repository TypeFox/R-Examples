psim.norm <-
function (q, mu, sig) {
	return(pnorm((q-mu)/sig/sqrt(mu^3*(1-mu)^3)))
}
