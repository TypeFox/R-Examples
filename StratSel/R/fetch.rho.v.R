fetch.rho.v <-
function(v,b){
	para <- b[length(b)]
	vari <- diag(v)[length(b)]
	br <- 2*(1/(1+exp(-para)))-1
	draws <- rnorm(1000,br,sqrt(vari))
	vv <- 2*(1/(1+exp(-draws)))-1
	varL <- var(vv)
	return(varL)
}
