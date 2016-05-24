loglikLOOCVcontourVAR1 <- function(lambdaAgrid, lambdaPgrid, Y, figure=TRUE, verbose=TRUE, ...){
	#############################################################################################
	# 
	# DESCRIPTION:
	# Evaluates the leave-one-out cross-validated log-likelihood of the VAR(1) model for a given grid of 
	# the ridge penalty parameters (lambdaA and lambdaO). The result is plotted as a contour plot, which 
	# facilitates the choice of optimal penalty parameters. The functions also works with a (possibly) 
	# unbalanced experimental set-up. The VAR(1)-process is assumed to have mean zero.
	#
	# ARGUMENTS:
	# -> lambdaAgrid   : Numeric of length larger than one. It contains the grid points corresponding to the lambdaA.
	# -> lambdaPgrid   : Numeric of length larger than one. It contains the grid points corresponding to the lambdaO.
	# -> Y             : Three-dimensional array containing the data. The first, second and third dimensions 
	#                    correspond to covariates, time and samples, respectively. The data are assumed to centered 
	#                    covariate-wise.
	# -> figure        : Logical, indicating whether the contour plot should be generated.
	# -> verbose       : Logical indicator: should intermediate output be printed on the screen?
	# -> ...           : Other arguments to be passed to loglikLOOCVVAR1.
	#
	# DEPENDENCIES:
	# require("rags2ridges")          # functions from package : loglikLOOCVVAR1
	#
	# NOTES:
	# ...
	#
	#############################################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(lambdaAgrid)) != "numeric"){ stop("Input (lambdaAgrid) is of wrong class.") }
	if (as.character(class(lambdaPgrid)) != "numeric"){ stop("Input (lambdaPgrid) is of wrong class.") }
	if (length(lambdaAgrid) < 2){ stop("Input (lambdaAgrid) is of wrong length.") }
	if (length(lambdaPgrid) < 2){ stop("Input (lambdaPgrid) is of wrong length.") }
	if (any(is.na(lambdaAgrid))){ stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") }
	if (any(is.na(lambdaPgrid))){ stop("Input (lambdaPgrid) is not a vector of non-negative numbers.") }
	if (any(lambdaAgrid <= 0)){ stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") }
	if (any(lambdaPgrid <= 0)){ stop("Input (lambdaPgrid) is not a vector of non-negative numbers.") }
	if (as.character(class(figure)) != "logical"){ stop("Input (figure) is of wrong class.") }
	if (as.character(class(verbose)) != "logical"){ stop("Input (verbose) is of wrong class.") }

	# evaluate cross-validated log-likelihood at all grid points
	lambdaAgrid <- sort(lambdaAgrid)
	lambdaPgrid <- sort(lambdaPgrid)

	llLOOCV <- matrix(NA, nrow=length(lambdaAgrid), ncol=length(lambdaPgrid))
	if (verbose){ cat("grid point:", "\n") }
	for (kA in 1:length(lambdaAgrid)){
		for (kP in 1:length(lambdaPgrid)){
			if (verbose){ 
		        cat(rep("\b", 100), sep ="")
				cat(paste("lambdaA=", lambdaAgrid[kA], "; lambdaP=", lambdaPgrid[kP], sep=""))
			}
			llLOOCV[kA, kP] <- loglikLOOCVVAR1(c(lambdaAgrid[kA], lambdaPgrid[kP]), Y, ...)
        	}     
    	}

	# plot contour
	if (figure){
		contour(lambdaAgrid, lambdaPgrid, -llLOOCV, xlab="lambdaA", ylab="lambdaP", main="cross-validated log-likelihood")
	}

	# return cross-validated 
	return(list(lambdaA=lambdaAgrid, lambdaP=lambdaPgrid, llLOOCV=-llLOOCV))
}



