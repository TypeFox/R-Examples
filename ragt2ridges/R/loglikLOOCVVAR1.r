loglikLOOCVVAR1 <- function(lambdas, Y, unbalanced=matrix(nrow=0, ncol=2), ...){
	#############################################################################################
	# 
	# DESCRIPTION:
	# Evaluation of the (minus) leave-one-out cross-validated log-likelihood of the VAR(1) model for 
	# given choices of the ridge penalty parameters (lambdaA and lambdaO). The functions also works 
	# with a (possibly) unbalanced experimental set-up. The VAR(1)-process is assumed to have mean zero.
	#
	# ARGUMENTS:
	# -> lambdas      : Numeric of length two. It contains the ridge penalty parameters to be used in the 
	#                   estimation of A and the Omega, the precision matrix of the errors (also called innovations). 
	# -> Y            : Three-dimensional array containing the data. The first, second and third dimensions 
	#                   correspond to covariates, time and samples, respectively. The data are assumed to 
	#                   centered covariate-wise. 
	# -> unbalanced   : A matrix with two columns, indicating the unbalances in the design. Each row represents 
	#                   a missing design point in the (time x individual)-layout. The first and second column 
	#                   indicate the time and individual (respectively) specifics of the missing design point.
	# -> ...          : Other arguments to be passed to ridgeVAR1.
	#
	# DEPENDENCIES:
	# require("rags2ridges")          # functions from package : ridgeVAR1
	#
	# NOTES:
	# ...
	#
	#############################################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(lambdas)) != "numeric"){ stop("Input (lambdas) is of wrong class.") }
	if (length(lambdas) != 2){ stop("Input (lambdas) is of wrong length.") }
	if (any(is.na(lambdas))){ stop("Input (lambdas) is not a vector of non-negative numbers.") }
	if (any(lambdas < 0)){ stop("Input (lambdas) is not a vector of non-negative numbers.") }
	if (as.character(class(unbalanced)) != "matrix"){ stop("Input (unbalanced) is of wrong class.") }    
	if (ncol(unbalanced) != 2){ stop("Wrong dimensions of the matrix unbalanced.") } 

	# determine leave-one-out scheme
	LOOscheme <- cbind(rep(2:dim(Y)[2], dim(Y)[3]), sort(rep(1:dim(Y)[3], dim(Y)[2]-1)))	
	if (nrow(unbalanced) > 0){
		LOO2unbalanced <- numeric()
		for (k in 1:nrow(unbalanced)){
			LOO2unbalanced <- c(LOO2unbalanced, which(apply(LOOscheme, 1, function(Y, Z){ all(Y == Z) }, Z=unbalanced[k,])))
		}
		LOOscheme <- LOOscheme[-LOO2unbalanced,]
	}

	loglik <- 0
	for (k in 1:nrow(LOOscheme)){
		# evaluate LOOCV estimates of VAR(1) parameters A and Se
		VAR1hat <- ridgeVAR1(Y, lambdas[1], lambdas[2], unbalanced=rbind(unbalanced, LOOscheme[k,,drop=FALSE]), ...)

		# evaluate LOOCV loglikelihood
		loglik <- loglik + .armaVAR1_loglik_LOOCVinternal(Y[,LOOscheme[k,1], LOOscheme[k,2]], Y[,LOOscheme[k,1]-1, LOOscheme[k,2]], VAR1hat$A, VAR1hat$P)
	}
	# return minus LOOCV loglikelihood
	return(-loglik)
}



