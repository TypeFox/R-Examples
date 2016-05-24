npSurface <-
function(fit){
	nt<-length(fit$fit)
	fit$n_regimes[1]+nt*(2+fit$n_regimes[2])
}
