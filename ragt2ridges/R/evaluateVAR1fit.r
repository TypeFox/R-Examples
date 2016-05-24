evaluateVAR1fit <- function(Y, A, SigmaE, unbalanced=NULL, diag=FALSE, fileType="eps", dir=getwd()){
	################################################################
	# 
	# DESCRIPTION:
	# ridge estimation of the parameters of the VAR(1) model
	#
	# ARGUMENTS: 
	# -> Y            : Three-dimensional array containing the data. The first, second and third 
	#                   dimensions correspond to covariates, time and samples, respectively. The 
	#                   data are assumed to centered covariate-wise. 
	# -> A            : Matrix A of regression parameters.
	# -> SigmaE       : Covariance matrix of the errors (innovations).
	# -> unbalanced   : A matrix with two columns, indicating the unbalances in the design. Each 
	#                   row represents a missing design point in the (time x individual)-layout. 
	#                   The first and second column indicate the time and individual (respectively) 
	#                   specifics of the missing design point.
	# -> diag         : Logical, should the diagonal be included in the evaluation of the fit. 
	# -> fileType     : Character specifying format in which figures should be save. Either 'pdf' or 'eps'.
	# -> dir          : Character specifying directory where plots should be saved.
	#
	# DEPENDENCIES:
	# library(marray)	    # functions: maPalette
	#
	# NOTES:
	# ...
	#
	################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(A)) != "matrix"){ stop("Input (A) is of wrong class.") }
	if (as.character(class(SigmaE)) != "matrix"){ stop("Input (SigmaE) is of wrong class.") }
	if (!isSymmetric(SigmaE)){ stop("Non-symmetric covariance matrix is provided.") }
	if (!all(eigen(SigmaE)$values > 0)){ stop("Non positive-definite covariance matrix is provided.") }
	if (nrow(A) != ncol(A)){ stop("Matrix A is not square.") }
	if (nrow(A) != nrow(SigmaE)){ stop("Dimensions covariance matrix and A do not match.") }
	if (!is.null(unbalanced) & as.character(class(unbalanced)) != "matrix"){ stop("Input (unbalanced) is of wrong class.") }    
	if (!is.null(unbalanced)){ if(ncol(unbalanced) != 2){ stop("Wrong dimensions of the matrix unbalanced.") } } 

	# extract dimension from data object
	nCovariates <- dim(Y)[1]
	nTimes <- dim(Y)[2]
	nSamples <- dim(Y)[3]

	# set profiles of missing (time, sample)-points to zero
	if (!is.null(unbalanced)){
		for (k in 1:nrow(unbalanced)){ Y[, unbalanced[k,1], unbalanced[k,2]] <- NA }
	}

	# obtain fits
	Yhat <- Y
	for (i in 1:nSamples){
		Yhat[,-1,i] <- A %*% Yhat[,-nTimes,i]
	}

	# evaluate fit A
	if (fileType=="pdf"){ pdf(file="plot_VAR1fit2observation.pdf") }
	if (fileType=="eps"){ setEPS(); postscript(file="plot_VAR1fit2observation.eps") }
	op <- par(pty="s")  
	ylims <- c(min(Y[,-1,], na.rm=TRUE), max(Y[,-1,], na.rm=TRUE)); xlims <- c(min(Yhat[,-1,], na.rm=TRUE), max(Yhat[,-1,], na.rm=TRUE))
	plot(y=Y, x=Yhat, col="white", xlim=xlims, ylim=ylims, ylab="observation", xlab="fit", main="assess fit (of A)")
	cols <- c("darkorange", "yellow", "orange", "red", "ivory", "firebrick", "indianred", "tomato", "khaki", "orangered", "chocolate")
	for (i in 1:nSamples){
		slhCols <- maPalette(low=paste(cols[i], "1", sep=""), high=paste(cols[(i %% 12) + 1], "3", sep=""), mid=NULL, k=nCovariates)
		for (j in 1:nCovariates){
			points(y=Y[j,-1,i], x=Yhat[j,-1,i], col=slhCols[j], pch=20)
		}
	}
	par(op)
	dev.off()

	# set profiles of missing (time, sample)-points to zero
	if (!is.null(unbalanced)){
		for (k in 1:nrow(unbalanced)){ Y[, unbalanced[k,1], unbalanced[k,2]] <- NA }
	}

	# obtain sample error covariance
	resSSY <- function(i, Z, Amat){ 
        	res <- cbind(Z[,1, i] - Amat %*% Z[,-dim(Z)[2], i])
	    	idContribute <- which(apply(!is.na(res), 2, all))
        	return(c(length(idContribute), as.numeric(res[,idContribute] %*% t(res[,idContribute]))))
	}    	 
	Se <- matrix(unlist(lapply(1:nSamples, resSSY, Z=Y, Amat=A)), nrow=1+nCovariates^2, byrow=FALSE)
   	Se <- matrix(rowSums(Se[-1,,drop=FALSE]), ncol=nCovariates) / sum(Se[1,,drop=FALSE])
   	
   	# evaluate fit of SigmaE
   	evaluateSfit(symm(solve(SigmaE)), Se, diag=diag, fileType=fileType, dir=dir)
}

