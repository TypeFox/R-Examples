`createGP` <-
function(X, Z, beta, a, meanReg, sig2, nugget, 
	param.names = 1:dim(X)[2], constantMean = 1) {

	l = list(beta,a)
	numDims = lapply(l, length)
	if (!all(unlist(lapply(numDims, "==", dim(X)[2])))) {
		stop("Length of  beta and/or a do not match dimensions of X")
	}
      
	gp = NULL
        gp$Z = Z
	gp$numObs = dim(X)[1]
	gp$numDim = dim(X)[2]
	gp$constantMean = constantMean

	if (constantMean == 1) gp$mu = matrix(meanReg, dim(X)[1])
	else gp$mu = cbind(  matrix(rep(1,dim(X)[1])), X) %*%meanReg

	gp$Bhat = meanReg
        gp$beta = beta
        gp$a = a
        gp$sig2 = sig2
	gp$params = param.names

        gp$X = X

        if (length(nugget) == 1) {
        	gp$invVarMatrix = try(solve(calcVarMatrix(X, beta, a, nugget,sig2,0, dim(X)[1])),TRUE)
		orig.nugget = nugget
		count = 1
		while(class(gp$invVarMatrix) == "try-error" && count < 1000) {
			nugget = 2*nugget + 1e-7
        		gp$invVarMatrix = try(solve(calcVarMatrix(X, beta, a, nugget,sig2,0, dim(X)[1])),TRUE)
			count = count + 1
		}
                if (count == 1000) {
			cat("\nERROR: final variance-covariance matrix is singular in R; try fitting setting min.nugget\n")
			return (NULL)
	        }
		if (count > 1) {
			cat("\nwarning: nugget has been changed to make var-cov matrix more stable\n")
			cat("from ")
			cat(orig.nugget)
			cat(" to ")
			cat(nugget)
			cat("\n")
		}	
		#nugget = orig.nugget
	}
        else {
                        gp$invVarMatrix = solve(calcVarMatrix(X, beta, a, 0,sig2,0, dim(X)[1])
                                + diag(as.vector(nugget), dim(X)[1]))
       }
        
	gp$nugget = nugget

	### gp$loglike = dmvnorm(t(Z), gp$mu, solve(gp$invVarMatrix), log=TRUE) ## numerical error
        gp$loglike = calcLogLikeManual(gp)
	gp$cv = NULL
	gp$cv <- CV(gp)
	attr(gp, "class") <- "gp"
        gp
}

