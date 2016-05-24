`summary.gp` <-
function(object, ...) {

	if (!is.null(object$numObs)) {
		cat("\nTotal observations = ")
		cat(object$numObs)
		cat(object$numReps)
		cat("\nDimensions = ")
		cat(object$numDim)
		cat("\n")
	}

	if (object$constantMean == 1) {
		cat("\nmu = ")
		cat(object$mu[1])
	}
	else {
		cat("\nmeanReg: ")
		cat(object$Bhat)
	}
	cat("\n")

	cat("sig2:\t")
	cat(object$sig2)
	cat("\n")

	if (!is.null(object$nugget)&length(object$nugget)==1) {
		cat("nugget:\t")
		cat(object$nugget)
		cat("\n")
	}


	cat("\nCorrelation parameters:\n\n")
	d = data.frame("beta" = object$beta, "a" = object$a)
	print(d)

	cat("\nLog likelihood = ")
	cat(object$loglike)
	cat("\n\n")

	if (!is.null(object$cv)) {
		RMSE = sqrt(mean((object$cv[,1] - object$Z)**2))
		RMaxSE = max((object$cv[,1] - object$Z)**2)	
		if (anyReps(object$X)) {
			v = varPerReps(object$X, object$Z)
			RMPE = sqrt(mean(v))
			cat ("CV RMSE (RMPE): ")
			cat (RMSE)
			cat( " (")
			cat(RMPE)
			cat(")\n")
			cat("CV RMaxSE (RMaxPE): ")
			index = which.max(v)
			cat (RMaxSE)
			cat (" (")
			cat(sqrt(v[index]))
			cat(")\n")
		}
		else {
			cat ("CV RMSE: ")
			cat (RMSE)
			cat("\n")

			cat("CV RMaxSE: ")
			cat (RMaxSE)
			cat("\n")
		}
	}
	

 
}

