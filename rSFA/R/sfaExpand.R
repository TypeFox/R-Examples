#TODO both functions should only be default functions, not hardcoded. Userdefined functions for expansion should be possible
# TODO: this has to be documented!

###################################################################################
#' Degree 2 Expansion
#'
#' Expand a signal in the space of polynomials of degree 2. 
#' This is the default expansion function used by rSFA. 
#'
#' @param sfaList 		A list that contains all information about the handled sfa-structure
#' @param DATA				Input data, each column a different variable
#'
#' @return expanded matrix \code{DATA}
#'
#' @seealso  \code{\link{sfa2}} \code{\link{nlExpand}} \code{\link{xpDim}}
#' @export
###################################################################################
sfaExpand <- function (sfaList, DATA){   
	n=customSize(DATA,2);
	m=customSize(DATA,1);
	DATA=cbind(DATA, matrix(0,m,sfaList$xpDimFun(n)-n))
	k=n+1;

	for (i in 1:n){
		len=n-i+1;
		#tmp=customRepmat(DATA[,i],1,len);  
		tmp=matrix(DATA[,i],m,len);#MZ, 11.11.12: speedfix
		DATA[,k:(k+len-1)]=tmp*DATA[,i:n];
		k=k+len;
	}	
	return(DATA)
}

###################################################################################
#' Degree 2 Dimension Calculation
#'
#' Compute the dimension of a vector expanded in the space of
#' polynomials of 2nd degree.
#'
#' @param n 			Dimension of input vector
#'
#' @return Dimension of expanded vector
#'
#' @seealso  \code{\link{sfa2}} \code{\link{sfaExpand}}
#' @export
###################################################################################
xpDim <- function (n){
	res=n+(n+1)*n/2;
	#res=choose(n*2,2)-1  #see paper by berkes: D= (N+d over d) -1
	#check if n can be removed once, and thus corresponding line in sfaExpand can also be removed
}

###################################################################################
#' Expand a signal in the for Nonlinear Expansion demo
#' 
#' This is an optional expansion function, calculating the expanded data in [x; x^4].
#'
#' @param sfaList 		A list that contains all information about the handled sfa-structure
#' @param DATA				Input data, each column a different variable
#'
#' @return expanded matrix \code{DATA}
#'
#' @seealso  \code{\link{sfa2}} \code{\link{sfaExpand}} \code{\link{nlDim}}
#' @export
#' @keywords internal
###################################################################################
nlExpand <- function (sfaList, DATA){    #example for user defined nlExpandFun (for expansion demo)
	return(cbind(DATA, DATA^4))
}

# % NL_DIM Compute the dimension of a vector expanded in space by 2.
# see nlExpand
###################################################################################
#' Custom Nonlinear Dimension Calculation
#'
#' Compute the dimension, used for the customexpansion demo.
#'
#' @param n 			Dimension of input vector
#'
#' @return Dimension of expanded vector
#'
#' @seealso  \code{\link{sfa2}} \code{\link{nlExpand}}
#' @export
#' @keywords internal
###################################################################################
nlDim <- function (n){ #example for user defined xpDimFun
	res=2*n
}