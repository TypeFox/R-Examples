`angledist` <-
function(distribution, distpars=NA){
	l <- list(distribution=distribution,chisq=NA, loglik=NA, AIC=NA,  distpars=distpars,
						fitmethod="Defined", fitdata = NA)
	class(l) <- "angledist"
	return(l)
}

