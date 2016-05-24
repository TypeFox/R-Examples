fpca <-
function(x, nbasisInit, propVar=.9, reconstruct=FALSE, varName=NULL, verbose=FALSE){

	interval <- 1:ncol(x)
	nbasis <- ifelse(missing(nbasisInit), ncol(x) / 4, nbasisInit)
	if(verbose) cat(nbasis, " Spline basis coefficients\n")
	bsp <- create.bspline.basis(c(1,ncol(x)), nbasis = nbasis)
	
	basis <- eval.basis(interval, bsp)
	fdObj <- Data2fd(argvals = interval, y=t(x), basisobj = bsp)




	fpca <- pca.fd(fdObj, nharm=nbasis, centerfns = TRUE)
	nrPC <- pmax(which(cumsum(fpca$varprop)>=propVar)[1], 2)
	if(verbose) cat(nrPC, "PCs selected\n")

	optimalDesign <- fpca$scores[,1:nrPC]

	str <- ifelse(is.null(varName), "PC", paste(varName, "PC", sep="_"))
	colnames(optimalDesign) <- paste(str, 1:nrPC, sep="")

	if(reconstruct){
		basisMean <- eval.basis(interval, fpca$meanfd$basis)
		meanFunction <- as.numeric(basisMean%*%fpca$meanfd$coefs)

		smoothData <- t(basis %*%  fpca$harmonics$coefs[,1:nrPC] %*% t(optimalDesign))
		smoothData <- t(apply(smoothData, MARGIN=1, FUN=function(z) z + meanFunction))

		lout <- list("design"=optimalDesign, "smoothData"=smoothData)
	}else{
		lout <- list("design"=optimalDesign, "smoothData"=NULL)
	}
	lout
}
