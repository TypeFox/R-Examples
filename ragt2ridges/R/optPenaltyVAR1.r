optPenaltyVAR1 <- function (Y, lambdaMin, lambdaMax, lambdaInit = (lambdaMin+lambdaMax)/2, 
                                  optimizer="nlm", ...){ 
	#####################################################################################################
	# 
	# DESCRIPTION: 
	# Ridge estimation of the parameters of the VAR(1) model. The log-likelihood is augmented with 
	# a ridge penalty for both parameters, A, the matrix of regression coefficients, and 
	# SigmaE, the inverse of the error variance. 
	# 
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The first, second and third dimensions correspond to 
	#                    covariates, time and samples, respectively. The data are assumed to centered covariate-wise.
	# -> lambdaMin     : Numeric of length two, containing the minimum values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter corresponding to the penalty on A, the matrix with regression 
	#                    coefficients, while the second parameter relates to the penalty on Omega, the precision matrix of the errors.
	# -> lambdaMax     : Numeric of length two, containing the maximum values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter corresponding to the penalty on A, the matrix with regression 
	#                    coefficients, while the second parameter relates to the penalty on Omega, the precision matrix of the errors.
	# -> lambdaInit    : Numeric of length two, containing the initial values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter corresponding to the penalty on A, the matrix with regression 
	#                    coefficients, while the second parameter relates to the penalty on Omega, the precision matrix of the errors.
	# -> optimizer     : A character : which optimization function should be used: "nlm" (default) or "optim"?
	# -> ...           : Additional arguments passed on to loglikLOOCVVAR1
	# 
	# DEPENDENCIES:
	# library(base)	        # functions: nlminb, constrOptim
	# library(ragt2ridges)  # functions: loglikLOOCVVAR1 and its dependencies.
	#
	# NOTES:
	# ....
	# 
	#####################################################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(lambdaMin)) != "numeric"){ stop("Input (lambdaMin) is of wrong class.") }
	if (length(lambdaMin) != 2){ stop("Input (lambdaMin) is of wrong length.") }
	if (any(is.na(lambdaMin))){ stop("Input (lambdaMin) does not comprise of positive numbers.") }
	if (any(lambdaMin < 0)){ stop("Input (lambdaMin) does not comprise of positive numbers.") }
	if (as.character(class(lambdaMax)) != "numeric"){ stop("Input (lambdaMax) is of wrong class.") }
	if (length(lambdaMax) != 2){ stop("Input (lambdaMax) is of wrong length.") }
	if (any(is.na(lambdaMax))){ stop("Input (lambdaMax) does not comprise of positive numbers.") }
	if (any(lambdaMax < 0)){ stop("Input (lambdaMax) does not comprise of positive numbers.") }
	if (any(lambdaMax <= lambdaMin)){ stop("Input (lambdaMax) must be larger (element-wise) than lambdaMin") }
	if (as.character(class(lambdaInit)) != "numeric"){ stop("Input (lambdaInit) is of wrong class.") }
	if (length(lambdaInit) != 2){ stop("Input (lambdaInit) is of wrong length.") }
	if (any(is.na(lambdaInit))){ stop("Input (lambdaInit) does not comprise of positive numbers.") }
	if (any(lambdaInit < 0)){ stop("Input (lambdaInit) does not comprise of positive numbers.") }
	if (any(lambdaInit <= lambdaMin)){ stop("Input (lambdaInit) must be larger (element-wise) than lambdaMin") } 
	if (any(lambdaInit >= lambdaMax)){ stop("Input (lambdaInit) must be smaller (element-wise) than lambdaMax") }

	# optimize LOOCV log-likelihood w.r.t. the penalty parameters
	if (optimizer=="optim"){    
		optLambdas <- constrOptim(lambdaInit, loglikLOOCVVAR1, grad=NULL, ui=rbind(diag(rep(1, 2)), diag(rep(-1, 2))), ci=c(lambdaMin, -lambdaMax), Y=Y, control=list(reltol=10^(-10)), ...)$par 
	}
	if (optimizer=="nlm"){    
		optLambdas <- nlminb(lambdaInit, loglikLOOCVVAR1, gradient=NULL, lower=lambdaMin, upper=lambdaMax, Y=Y, control=list(rel.tol=10^(-10)), ...)$par 
	}
	return(optLambdas)
}


