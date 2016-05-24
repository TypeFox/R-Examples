#'
#' The C++ code, of this package take to run from the data two main information.
#' the coefficient in the basis expansion of the functional data, and the inner product between theses basis.
#' cppUnidata made this task in the univariate case.
#' 
#' @param fd the functional data 
#' 
#' @return fdData a list that containing the transpose of the coefficients matrix and the inner 
#' product between these basis functions. 
#' 
cppUniData <- function(fd) {
	coefs=t(fd$coefs)
	basisProd=inprod(fd$basis,fd$basis)
	fdData=list(coefs=coefs,basisProd=basisProd)
	return(fdData)
}


#'
#' The C++ code, of this package take to run from the data two main information.
#' the coefficient in the basis expansion of the functional data, and the inner product between theses basis.
#' cppMultidata made this task in the multivariate case.
#' 
#' @param mfd a list containing all dimension of the data
#' 
#' @return fdData a list that containing the concatenation of the transpose of the coefficients matrix 
#' and a block diagonal matrix where each block is the inner product between the basis functions.
#' 
cppMultiData <-function(mfd) {
	#mfd is a list of functional data (mfd=list(fd_1,fd_2,...,fd_dim))
	nbasis=c();
	dim=length(mfd)
	for (i in 1:dim) {
		nbasis[i]=mfd[[i]]$basis$nbasis;
	}
	basisProd=matrix(0,nrow=sum(nbasis),ncol=sum(nbasis))
	nobs=nrow(t(mfd[[1]]$coefs))
	coefs=matrix(0,nrow=nobs,ncol=sum(nbasis))
	i=0;
	iter=1;
	while (iter <= dim) {
		a=inprod(mfd[[iter]]$basis,mfd[[iter]]$basis);
		n=nbasis[iter];
		basisProd[(i+1):(i+n),(i+1):(i+n)]=a;
		coefs[,(i+1):(i+n)]=t(mfd[[iter]]$coefs)
		i=i+n;
		iter=iter+1;
	}
	fdData=list(coefs=coefs,basisProd=basisProd)
	return(fdData)
}
