summary.blca.gibbs <-
function(object, ...){
	sum1<- c(nrow(object$samples$classprob), object$burn.in, object$thin, mean(object$samples$logpost), object$logpost, object$DIC, object$BICM, object$AICM)
	names(sum1)<- c("IterNumber", "Burn-in", "ThiningRate", "Dbar", "Dhat", "DIC", "BICM", "AICM")
	
	object$method<- "Gibbs Sampling"
	object$printnames<- c("Chain Lengths:", "Burn-In:", "Thinning Rate:", "Log-Posterior Mean:", "Log-Posterior Mode:", "DIC:", "BICM:", "AICM:")
	object$sum1<- sum1

	NextMethod("summary")
	}
