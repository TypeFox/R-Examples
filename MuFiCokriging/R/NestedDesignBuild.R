NestedDesignBuild <- function(design = NULL){
	
	if(identical(design,NULL)){stop("The user must enter nlevel design sets")}

	nlevel <- length(design)

	PX <- list()
	indices <- list()

	PX[[1]] <- design[[1]]
	nPX <- dim(PX[[1]])[1]
	dPX <- dim(PX[[1]])[2]
	for(i in 2:nlevel){
		PX[[i]] <- design[[i]]
		nPX <- c(nPX , dim(PX[[1]])[1])
		dPX <- c(dPX , dim(PX[[1]])[2])
		if(nPX[i] > nPX[i-1]){stop("The number of experiments at level i must be lower than the one at level i-1")}
		if(dPX[i] > dPX[i-1]){stop("The number of dimensions is not consistent")}
	}

	for(i in (nlevel-1):1){
		SB <- SubstDesign(PX[[i+1]],PX[[i]])
		PX[[i]] <- SB$PX
		n <- dim(SB$PX)[1]
		indices[[i]] <- seq(n-SB$le+1,n,by=1)
	}

	PX <- NestedDesign(PX[[1]], nlevel = nlevel , indices = indices)
	return(PX)
}