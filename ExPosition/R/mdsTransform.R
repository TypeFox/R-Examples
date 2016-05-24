mdsTransform <- function(D,masses){
	DATA_dimensions <- dim(D)
	#do this every time.
	Mrepmat <- matrix(masses,nrow=nrow(D),ncol=ncol(D))
	if(is.null(dim(masses))){ # it is a vector; with new way, it is always a vector
		BigXi <- diag(DATA_dimensions[1]) - (matrix(1,DATA_dimensions[1],1) %*% masses)
	}else{#ths forces a matrix to be a vector this needs to be better.
		BigXi <- diag(DATA_dimensions[1]) - (matrix(1,DATA_dimensions[1],1) %*% diag(masses))
	}
	S <- -.5 * sqrt(Mrepmat) * BigXi %*% D %*% t(BigXi) * sqrt(t(Mrepmat))
	rownames(S) <- rownames(D)
	colnames(S) <- colnames(D)
	return(S)
}