#'
#' This function is used in mfpca to cut the harmonic coefficient matrix. In fact
#' mfpca call a c++ mfpca which return (among others) the matrix of the coefficients of the harmonics 
#' in the basis expantion. But in the multivariate case the coefficient matrix for all dimension
#' are store in the same matrix, so we use this function to store each dimension in the specific matrix.     
#' 
#' @title Separates the matrices of the coefficients of harmonics
#' 
#' @param harms the coefficients matrix, see the description for more details
#' @param nbasis vector containing the number of basis for each dimension
#' 
#' @return a list of length dimension of the data which contain the harmonic coefficients matrix
#' for each dimension
#' 
harmsCut <-function(harms,nbasis){
	dimFd=length(nbasis)
	i=0
	iter=1
	harmsCut=list()
	while (iter <= dimFd){
		start=i+1;
		end=i+nbasis[[iter]]
	 	harmsCut[[iter]]=harms[start:end,]
	 	i=end
	 	iter=iter+1
	}
	return(harmsCut)
}


