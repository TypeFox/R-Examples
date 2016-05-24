loglikVAR1 <- function(Y, A, P, unbalanced=matrix(nrow=0, ncol=2)){
	##############################################################################
	# 
	# DESCRIPTION:
	# Log-likelihood of the VAR(1) model specified by the supplied parameters.
	#
	# ARGUMENTS: 
	# -> Y            : Three-dimensional array containing the data. The first, second and third 
	#                   dimensions correspond to covariates, time and samples, respectively. The 
	#                   data are assumed to centered covariate-wise. 
	# -> A            : Matrix A of regression parameters.
	# -> P            : Inverse error covariance matrix.
	# -> unbalanced   : A matrix with two columns, indicating the unbalances in the design. Each 
	#                   row represents a missing design point in the (time x individual)-layout. 
	#                   The first and second column indicate the time and individual (respectively) 
	#                   specifics of the missing design point.
	#
	# DEPENDENCIES:
	# ...
	#
	# NOTES:
	# ...
	#
	##############################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(A)) != "matrix"){ stop("Input (A) is of wrong class.") }
	if (as.character(class(P)) != "matrix"){ stop("Input (P) is of wrong class.") }
	# if (!isSymmetric(P)){ stop("Non-symmetric precision matrix is provided.") }
	if (!all(eigen(P)$values > 0)){ stop("Non positive-definite precision matrix is provided.") }
	if (nrow(A) != ncol(A)){ stop("Matrix A is not square.") }
	if (nrow(A) != nrow(P)){ stop("Dimensions precision matrix and A do not match.") }
	if (as.character(class(unbalanced)) != "matrix"){ stop("Input (unbalanced) is of wrong class.") }    
	if (ncol(unbalanced) != 2){ stop("Wrong dimensions of the matrix unbalanced.") } 

	# set profiles of missing (time, sample)-points to missing
	Y <- .armaVAR_array2cube_withMissing(Y, unbalanced[,1], unbalanced[,2])

	# obtain loglikelihood contribution of residuals 
	LL <- .armaVAR1_loglik(Y, A, P)
	return(LL)
}

