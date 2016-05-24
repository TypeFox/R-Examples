`fitalldistributions` <-
function(angles,
								fitmethod="loglik", 
								distributions = c("twoparbeta", "ellipsoid", "rotatedell", 
					                              "planophile","erectophile","plagiophile",
					                              "extremophile","spherical","uniform"),
								...){

	fitobjs <- list()
	chisqs <- AICs <- c()
	for(i in 1:length(distributions)){
		d <- distributions[i]
		fitobjs[[i]] <- fitdistribution(angles, distribution=d, fitmethod=fitmethod, ...)
		
		if(fitmethod=="chisq")chisqs[i] <- fitobjs[[i]]$chisq
		if(fitmethod=="loglik")AICs[i] <- fitobjs[[i]]$AIC

	}
	
	names(fitobjs) <- distributions
	bestfit <- rep("", length(distributions))
	
	if(fitmethod=="chisq")bestfit[which.min(chisqs)] <- "best fit"
	if(fitmethod=="loglik")bestfit[which.min(AICs)] <- "best fit"
	
if(fitmethod=="chisq"){

	ans <- list(allfits = fitobjs,
				fitsummary = data.frame(distribution=distributions, chisq=chisqs, bestfit=bestfit))
	class(ans) <- "angledistlist"
	return(ans)
}

if(fitmethod=="loglik"){
	ans <- list(allfits = fitobjs,
				fitsummary = data.frame(distribution=distributions, AIC=AICs, bestfit=bestfit))
	class(ans) <- "angledistlist"
	return(ans)
}

}

