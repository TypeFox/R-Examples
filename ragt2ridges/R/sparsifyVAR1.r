sparsifyVAR1 <- function(A, SigmaE, threshold=c("absValue", "localFDR", "top"), absValueCut=0.25, FDRcut=0.8, top=10, statistics=FALSE, verbose=FALSE){
	#####################################################################################################
	#
	# DESCRIPTION: 
	# Function that sparsifies/determines support of the matrix A of a VAR(1) model
	# Support can be determined by absolute value thresholding or by local FDRs thresholding
	# One can also choose to threshold based on the top X of absolute value of the elements of A
	# Function is to some extent a wrapper around certain 'fdrtool' functions
	#
	# ARGUMENTS:
	# -> A           : (Possibly shrunken) matrix with regression coefficients
	# -> SigmaE      : (Possibly shrunken) error covariance matrix
	# -> threshold   : Signifies type of thresholding
	# -> absValueCut : Cut-off for elements selection based on absolute value. 
	#                  thresholding (of A). Only when threshold = 'absValue'. Default = .25
	# -> FDRcut      : Cut-off for element selection based on local FDR thresholding on test statistics.
	#                  Only when threshold = 'localFDR'. Default = .9
	# -> top         : element selection based on retainment 'top' number of absolute values of A.
	#                  Only when threshold = 'top'. Default = 10
	# -> statistics  : logical indicating if test statistics should be returned
	#                  Only when threshold = 'localFDR'. Default = FALSE
	# -> verbose     : Logical indicating if intermediate output should be printed on screen
	#                  Only when threshold = 'localFDR'. Default = TRUE
	#
	# DEPENDENCIES:
	# require("fdrtool")          # functions from package : fdrtool
	#
	# NOTES:
	# For the sparsification of the (inverse of the) SigmaE matrix, confer 'sparsifyS'
	# 
	#####################################################################################################

	# input check
	if (as.character(class(A)) != "matrix"){ stop("Input (A) is of wrong class.") }
	if (nrow(A) != ncol(A)){ stop("Matrix A is not square.") }	
	if (as.character(class(threshold)) != "character"){ stop("Input (threshold) is of wrong class.") }
	if (missing(threshold)){ stop("Need to specify type of sparsification ('absValue' or 'localFDR' or 'top')") }
	if (!(threshold %in% c("absValue", "localFDR", "top"))){ stop("Input (threshold) should be one of {'absValue', 'localFDR', 'top'}") }
	if (as.character(class(verbose)) != "logical"){ stop("Input (verbose) is of wrong class.") }

	# top elements of A
	if (threshold == "top"){
		# extra input checks for this option
		if (class(top) != "numeric"){ stop("Input (top) is of wrong class") } 
		if (length(top) != 1){ stop("Input (top) must be a scalar") }
		if (!.is.int(top)){ stop("Input (top) should be a numeric integer") }
		if (top <= 0){ stop("Input (top) must be strictly positive") } 
		if (top >= (ncol(A))^2){ stop("Input (top) must be smaller than the number of elements of the input matrix A") }
		
		# determine threshold for top
		absValueCut <- sort(abs(A), decreasing = TRUE)[ceiling(top)]
		threshold   <- "absValue"
	}

	# absolute value criterion
	if (threshold == "absValue"){
		# extra input checks for this option
		if (class(absValueCut) != "numeric"){ stop("Input (absValueCut) is of wrong class") } 
		if (length(absValueCut) != 1){ stop("Input (absValueCut) must be a scalar") }  
		if (absValueCut <= 0){ stop("Input (absValueCut) must be positive") }
        
		# selection
		nonzeros <- which(abs(A) >= absValueCut, arr.ind=TRUE)
		zeros <- which(abs(A) < absValueCut, arr.ind=TRUE)
		H0statA <- NULL
	}

	# local FDR values 
	if (threshold=="localFDR"){
		# extra input checks for this option
		if (as.character(class(SigmaE)) != "matrix"){ stop("Input (SigmaE) is of wrong class.") }
		if (!isSymmetric(SigmaE)){ stop("Non-symmetric covariance matrix is provided.") }
		if (!all(eigen(SigmaE, only.values=TRUE)$values > 0)){ stop("Non positive-definite covariance matrix is provided.") }
		if (nrow(A) != nrow(SigmaE)){ stop("Dimensions covariance matrix and A do not match.") }
		if (as.character(class(FDRcut)) != "numeric"){ stop("Input (FDRcut) is of wrong class.") }
		if (length(FDRcut) != 1){ stop("Input (FDRcut) is of wrong length.") }
		if (is.na(FDRcut)){ stop("Input (FDRcut) is not a positive number.") }
		if (FDRcut <= 0){ stop("Input (FDRcut) is not a positive number.") }
		if (FDRcut >= 1){ stop("Input (FDRcut) should be smaller than one.") }
		if (as.character(class(statistics)) != "logical"){ stop("Input (testStat) is of wrong class.") }
	
		# calculate the variance of Y	
		Syy <- SigmaE
		for (tau in 1:1000) {
			Atau <- A %^% tau
			Syy <- Syy + Atau %*% SigmaE %*% t(Atau)
			if (max(abs(Atau)) < 10^(-20)){ break }
		}

		# test statistics for elements of A: beta-hat / s.e.(beta-hat)
		H0statA <- as.numeric(A / sqrt(t(outer(diag(solve(Syy)), diag(SigmaE)))))
    
		# 'testing'
		lFDRs <- 1 - fdrtool(H0statA, "normal", cutoff.method="locfdr", plot=verbose, verbose=verbose)$lfdr
		nonzeros <- which(matrix(lFDRs, ncol=ncol(A), byrow=FALSE) > FDRcut, arr.ind=TRUE)
		zeros <- which(matrix(lFDRs, ncol=ncol(A), byrow=FALSE) <= FDRcut, arr.ind=TRUE)
	}    
	
	if (verbose){
        cat("-> Retained elements: ", nrow(nonzeros), "\n")
        cat("-> Corresponding to", round(nrow(nonzeros)/(nrow(nonzeros) + nrow(zeros)), 4) * 100, "% of possible elements \n")
    }
    
	# return the elements of 	
	if (!statistics){ H0statA <- NULL }
	return(list(nonzeros=nonzeros, zeros=zeros, statistics=H0statA))	    
}


