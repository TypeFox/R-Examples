mdsSupplementalElementsPreProcessing <- function(SUP.DATA=NULL,D=NULL,M=NULL){

	if(is.null(SUP.DATA) || is.null(D)){
		stop('You must provide supplemental and active data.')
	}
	if(nrow(D)!=ncol(D)){
		stop('D dims do not match.')
	}
	if(sum(diag(D))!=0){
		stop('D diag is not 0.')
	}	
	if(nrow(SUP.DATA)!=ncol(SUP.DATA)){
		stop('SUP.DATA dims do not match.')
	}
	if(sum(diag(SUP.DATA))!=0){
		stop('SUP.DATA diag is not 0.')
	}
	if(nrow(SUP.DATA)!=nrow(D)){
		stop('SUP.DATA and D dims do not match.')
	}
	if(is.null(M)){
		M <- rep(1/nrow(D),nrow(D))
	}

	#M needs to be a vector here. D needs to be the original distance matrix. SUP.DATA needs to be supp. data in a dist matrix.
	BigXi <- diag(nrow(D)) - (matrix(1,nrow(D),1) %*% M)
	Mrepmat <- matrix(M,nrow=nrow(D),ncol=ncol(D))
	sup.subtract <- (SUP.DATA - (D %*% as.matrix(M) %*% matrix(1,1,length(M))))
	
	return(-0.5 * sqrt(Mrepmat) * BigXi %*% sup.subtract %*% t(BigXi) * sqrt(t(Mrepmat)))
}