fetch.rho.b <-
function(b){
	para <- b[length(b)]
	b[length(b)] <- 2*(1/(1+exp(-para)))-1
	return(b)
}
