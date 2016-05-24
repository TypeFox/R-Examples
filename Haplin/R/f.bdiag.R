f.bdiag <- function(matlist){
##
##
.d <- sapply(matlist, function(x)dim(x)[1])
.dnames1 <- unlist(lapply(matlist, function(x)dimnames(x)[[1]]))
.dnames2 <- unlist(lapply(matlist, function(x)dimnames(x)[[2]]))

.mat <- matrix(0, ncol = sum(.d), nrow = sum(.d), dimnames = list(.dnames1, .dnames2))




.f.start.stop <- function(d){
## COMPUTE STARTPOINTS AND ENDPOINTS IN AN EXTRACTION SEQUENCE
## NO, SIMPLER JUST TO MAKE THE ACTUAL SEQUENCE:
	.dc1 <- cumsum(d) # END POINTS
	.dc0 <- c(0, .dc1[-length(.dc1)]) + 1 # START POINTS
	.ut <- NULL
	for (i in seq(along = .dc0)){
		.ut[[i]] <- .dc0[i]:.dc1[i]
	}
	###return(data.frame(start = .dc0, stop = .dc1))
	return(.ut)
}


.f.insert <- function(matfull, matsublist){
## FUNCTION FOR INSERTING A LIST OF SUBMATRICES (BLOCK DIAGONALS) INTO A LARGER MATRIX
	.d <- sapply(matsublist, function(x)dim(x)[1])
	if(dim(matfull)[1] != sum(.d))stop()
	#
	.st <- .f.start.stop(.d)
	for(i in seq(along = matsublist)){
		matfull[.st[[i]], .st[[i]]] <- matsublist[[i]]
	}
	return(matfull)
}


.f.extract <- function(matfull, d){
## FUNCTION FOR EXTRACTING A BLOCK-DIAGONAL FROM A LARGER MATRIX
	if(sum(d) != dim(matfull)[1]) stop()
	.ut <- vector(length(d), mode = "list")
	.st <- .f.start.stop(d)
	for (i in seq(along = d)){
		.ut[[i]] <- matfull[.st[[i]], .st[[i]], drop = F]
	}
	return(.ut)
}

.mat <- .f.insert(.mat, matlist)
return(.mat)


}

