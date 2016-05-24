estimable_R <- function(design, v, model, C, verbose=0) {
  model <- getModelNr(model)
	if(missing(C)) {
		Csub <- contrMat(n=rep(1, v), type="Tukey")
		class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
		C <- appendZeroColumns(Csub, model, v)
	}
	rcDesign <- rcd(design, v=v, model=model)
	Xr <- rcdMatrix(rcDesign, v, model)
	H <- linkMatrix(model, v)
	X <- Xr %*% H
	Z <- getZ(s=dim(design)[2],p=dim(design)[1])
	X <- cbind(X, Z)  
	XX <- t(X) %*% X
	C2 <- cbind(C, matrix(0, dim(C)[1], dim(Z)[2]))
	if (verbose) {            
		print(rcDesign)
		print(Xr)
		print(abs(C2 %*% ginv(XX) %*% XX-C2))
	}
	return(isTRUE(all.equal(C2 %*% ginv(XX) %*% XX, C2, check.attributes=FALSE, check.names=FALSE))) # Criterion to test whether - see Theorem \ref{thr:estimable} of vignette.
}

estimable <- function(design, v, model, C, verbose=0) {
  model <- getModelNr(model)
	if(missing(C)) {
		Csub <- contrMat(n=rep(1, v), type="Tukey")
		class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
		C <- appendZeroColumns(Csub, model, v)
	}
	rcDesign <- rcd(design, v=v, model=model)
	linkM <- linkMatrix(model, v)
	Z <- getZ(s=dim(design)[2],p=dim(design)[1])
	return(.Call( "estimable2R", rcDesign, v, model, linkM, C, Z, verbose, PACKAGE = "Crossover" ))    
}

meanEff <- function(eff) {
  mean(eff[upper.tri(eff)])
}