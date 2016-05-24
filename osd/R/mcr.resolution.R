mcr.resolution <-
function(D, k)
{

	S <- matrix(0, ncol=k, nrow=ncol(D))
	if(is.matrix(S)==F) dim(S) <- c(length(S),1)

	#PCA initialization
	pca.ini <- prcomp(D, center=F)

	w.sign <- diag(sign(as.matrix(pca.ini$x)[apply(as.matrix(abs(pca.ini$x)),2, which.max),]))
	w.sign[w.sign==0] <- -1
	pca.ini$x <- sweep(as.matrix(pca.ini$x),2,w.sign,"*")
	pca.ini$x[pca.ini$x<0] <- 0
	
	pca.ini.maxs <- apply(pca.ini$x,2,which.max)
	k.inds <- which(duplicated(pca.ini.maxs)==F)[1:k]

	C <- pca.ini$x[,k.inds]
	
	linear.transformation.gain <- max(D, na.rm=T)
	D <- D/linear.transformation.gain
	C <- C/max(C, na.rm=T)

	MCR.res <- mcrnonneg(D,C,S)
	MCR.res$resC <- MCR.res$resC*linear.transformation.gain
			
	return(list(C=MCR.res$resC,S=MCR.res$resS))
}
