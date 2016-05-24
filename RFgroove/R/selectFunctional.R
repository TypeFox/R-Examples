selectFunctional <-
function(FDlist, ydata, normalize=TRUE, dimensionReductionMethod=c("fpca", "wave"), nbasisInit, verbose=TRUE, ...){
	
	n <- nrow(FDlist[[1]]) # nr of examples
	N <- ncol(FDlist[[1]]) # dimension of each curve
	p <- length(FDlist) # nr of functional variables
	J <- log2(N) # nr of wavelet levels
	varNames <- names(FDlist) # names of the functional covariates

	if(dimensionReductionMethod=="wave" & J != round(J))
		stop("N must be a power of two for wavelet dimension reduction.")

	if(dimensionReductionMethod=="fpca" & missing(nbasisInit))
		stop("For 'fpca' reduction dimension method, the parameter 'nbasisInit' must be given.")


	# 1- Normalization + Projection : In = FDlist, Out = (designMatrix, nvarGroup)
	if(dimensionReductionMethod=="wave"){
		reductionOneFd <- function(j, ...){
			FDmatrix <- FDlist[[j]]
			varName <- varNames[[j]]
			if(verbose) cat(varName,"\n")

			if(normalize){
				cat("Normalization\n")
				meanNorm <- mean(apply(FDmatrix, MARGIN=1, FUN=function(v) sqrt(v%*%v)))
				FDmatrix <- t(apply(FDmatrix, MARGIN=1, FUN=function(v) v / meanNorm))
			}
			hardThresholding(xdata=FDmatrix, verbose=verbose, varName=varName, ...)$estimatedDesign
		}
		designMatrixList <- sapply(1:p, FUN=reductionOneFd)
	}
	if(dimensionReductionMethod=="fpca"){
		reductionOneFd <- function(j, ...){
			FDmatrix <- FDlist[[j]]
			varName <- varNames[[j]]
			if(verbose) cat(varName,"\n")

			if(normalize){
				if(verbose) cat("Normalization\n")
				meanNorm <- mean(apply(FDmatrix, MARGIN=1, FUN=function(v) sqrt(v%*%v)))
				FDmatrix <- t(apply(FDmatrix, MARGIN=1, FUN=function(v) v / meanNorm))
			}
			fpca(x=FDmatrix, nbasisInit=nbasisInit, verbose=verbose, varName=varName, ...)$design
		}
		designMatrixList <- sapply(1:p, FUN=reductionOneFd)
	}

	designMatrix <- nvarGroup <- numeric(0)
	for(j in 1:length(designMatrixList)){
		designMatrix <- cbind(designMatrix, designMatrixList[[j]])
		nvarGroup <- c(nvarGroup, ncol(designMatrixList[[j]]))
	}
	names(nvarGroup) <- varNames


	# 2- Selection
	selectGroup(designMatrix, ydata, varNames, nvarGroup, verbose=verbose, ...)
}