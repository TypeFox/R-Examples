mlegp2 <- function(ZZ, XX, orig.XX, nugget, nugget.known, reps, un.sum, numEstimates,constantMean, simplex.ntries, simplex.maxiter, simplex.reltol, BFGS.maxiter, BFGS.tol, BFGS.h, min.nugget, verbose, param.names = param.names) {

	    Z = ZZ$Z
	    i = ZZ$i	
	    seed = ZZ$seed

	    if (verbose > 0) {
		    cat("========== FITTING GP # "); cat(i);
		    cat(" ==============================\n")
	    }
		
	   if (nugget.known == 1 ) { 
		if(is.null(nugget)) {    ## calculate the BLUE of the nugget; each sample variance scaled by (n_i - 1) / (N-k)
			## check that reps > 1 for at least one
			if (max(reps) < 2) {
				stop("error: cannot estimate nugget when replicate runs are not available; set nugget to its known value or include replicate observations")
			}
			index = reps > 1
			nugget = sum( un.sum$uniqueVar[index,i]*(reps[index]-1) / (sum(reps[index]) - sum(index)))
		}
		if (length(nugget) > 1) {   ## nugget matrix is for unique obs
			nugget = uniqueSummary(orig.XX, matrix(nugget))$uniqueMeans	
		}
		nugget = nugget / reps 
	}
	if (is.null(nugget) || (length(nugget) == 1 && nugget != 0)) {
		if (anyReps(XX)) nugget = estimateNugget(XX,Z )	
		
	}

	    simplex.abstol = -99999

	    success = 0
	    estimates = rep(0,numEstimates)
 	    returnFromC = .C("fitGPfromR", as.double(XX), as.integer(nrow(XX)), as.integer(ncol(XX)),
		as.double(Z), as.integer(nrow(Z)), 
		as.integer(constantMean), 
		as.integer(simplex.ntries), as.integer(simplex.maxiter), 
			as.double(simplex.abstol), as.double(simplex.reltol),
		as.integer(BFGS.maxiter), as.double(BFGS.tol), as.double(BFGS.h), 
		as.integer(seed), as.double(nugget), as.integer(length(nugget)), as.double(min.nugget),
		estimates = as.double(estimates), verbose = as.integer(verbose), nugget.known = as.integer(nugget.known), success = as.integer(success), 
		PACKAGE="mlegp")

	    if (returnFromC$success != 0) {
		cat("ERROR: GP cannot be created\n")
	        return (NULL)
            }	
	    estimates = returnFromC$estimates

    	    numParams = dim(XX)[2]
	    regSize = 1
	    meanReg = estimates[1]
	    if (constantMean == 0) {
		regSize = dim(XX)[2] + 1
		meanReg = estimates[1:regSize]
	    }

 	    ## remove meanReg params; now we have correlation params, sig2, and (possibly) nugget
	    estimates = estimates[(regSize+1):length(estimates)]   ## remove mean reg params

	    beta = estimates[1:numParams]
	    a = rep(2, numParams)
	    sig2 = estimates[numParams+1]

	    if (is.null(nugget)) nugget = 0

	    if (length(nugget) > 1 && nugget.known == 0) {
		if (verbose > 0) {
			cat(paste("nugget matrix will be scaled by: ", estimates[numParams+2]))
			cat("\n")
		}
		nugget = nugget * estimates[numParams+2]
	    }
	    else {
		if (nugget > 0 && nugget.known == 0) {
			nugget = estimates[numParams+2]
		}
	    }
	

 	    if (is.null(param.names)) param.names = paste("p",1:dim(XX)[2],sep="")

	   nugget = nugget + min.nugget
	
	    if (verbose > 0) cat("creating gp object...")
	    fit1 = createGP(XX, as.matrix(Z), beta, a, meanReg, sig2, nugget, 
		param.names = param.names, constantMean = constantMean)
	    if(is.null(fit1)) {
		cat("GP #"); cat(i); cat(" cannot be created...exiting mlegp\n")
		return (NULL)
	    } 
		
	   ## report nugget for single observation
	   if (nugget.known ==1 ) {
		fit1$nugget = fit1$nugget*reps  
		uniqueNug = unique(fit1$nugget)
		if(length(uniqueNug) == 1) fit1$nugget = as.numeric(uniqueNug)  
	    }	
	    if (verbose > 0) cat("...done\n\n")

	    return(fit1)
}






