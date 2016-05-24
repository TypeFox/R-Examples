summary.MuFicokm <- function(object, CrossValidation = FALSE, ...){
	response <- object$response
	Nestdesign <- object$Dnest
	nlevel <- object$nlevel
	
	n <- dim(as.matrix(ExtractNestDesign(Nestdesign,nlevel)))[1]
	if(CrossValidation){
		errLOO <- apply(matrix(1:Nestdesign$n),1,function(x) CrossValidationMuFicokmAll(object,x)$CVerrall)
		varLOO <- apply(matrix(1:Nestdesign$n),1,function(x) CrossValidationMuFicokmAll(object,x)$CVvarall)

		RMSE.LOO <- mean(errLOO^2)
		Q2.LOO <- 1 - sum(errLOO^2)/sum((response[[nlevel]]-errLOO-mean(response[[nlevel]]))^2)
		Std.RMSE.LOO <- mean(errLOO^2/varLOO)
	}

	name.cov <- list()
	cov.estimate <- list()
	sigma2.estimate <- list()
	rho.estimate <- list()
	beta.estimate <- list()

	objecti <- object$cok[[1]]
	name.cov[[1]] <- objecti@covariance@name
	cov.estimate[[1]] <- covparam2vect(objecti@covariance)
	sigma2.estimate[[1]] <- objecti@covariance@sd2
	beta.estimate[[1]] <- objecti@trend.coef

	for( i in 2:nlevel){
		objecti <- object$cok[[i]]
		name.cov[[i]] <- objecti@covariance@name
		cov.estimate[[i]] <- covparam2vect(objecti@covariance)
		sigma2.estimate[[i]] <- objecti@covariance@sd2
		mu.estimate <- objecti@trend.coef
		db <- dim(objecti@AR.F)[2]
		rho.estimate[[i-1]] <- mu.estimate[1:db]
		beta.estimate[[i]] <- mu.estimate[(db+1):length(mu.estimate)]
	}
	
	
	## print
	cat("Level 1: parameter estimation \n","\n",
	    	"Covariance type:",name.cov[[1]],"\n",
		"Variance estimation:",sigma2.estimate[[1]],"\n",
		"Trend estimation:",beta.estimate[[1]],"\n",
		"Correlation length estimation:",cov.estimate[[1]],"\n","\n",sep=" ")
		## print
	for(i in 2:nlevel){
	cat("Level ",i,": parameter estimation \n","\n",
	    	"Covariance type:",name.cov[[i]],"\n",
		"Variance estimation:",sigma2.estimate[[i]],"\n",
		"Trend estimation:",beta.estimate[[i]],"\n",
		"Adjustment estimation:",rho.estimate[[i-1]],"\n",
		"Correlation length estimation:",cov.estimate[[i]],"\n","\n",sep=" ")
	}
	if(CrossValidation){
		cat("Leave One Cross Validation \n","\n",
		    	"RMSE:",RMSE.LOO,"\n",
			"Std RMSE:",Std.RMSE.LOO,"\n",
			"Q2:",Q2.LOO,"\n","\n",sep="")
	}
	return(list(
	CovNames = name.cov,
	Cov.Val = cov.estimate ,
	Var.Val = sigma2.estimate ,
	Rho.Val = rho.estimate ,
	Trend.Val = beta.estimate
	))
}





