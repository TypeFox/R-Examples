###################################################################################
#' Create a new covariance object. 
#' 
#' @param dim	dimension of the covariance object list to be created
#'
#' @return list containing the COVMTX, avg and tlen
#'
#' @seealso  \code{\link{lcovUpdate}}
#' @export
#' @keywords internal
###################################################################################
lcovCreate <- function (dim){
	lcov=list()
	lcov$COVMTX=matrix(0,dim,dim);
	lcov$avg=matrix(0,1,dim);
	lcov$tlen=0;
  	return(lcov)
}

###################################################################################
#' Update a covariance object 
#' 
#'  Updates the covariance object
#'  referenced by lcov with a new chunk of data DATA. 
#'
#' @param lcov 			A list that contains all information about the handled covariance-structure
#' @param DATA			Input Data: must be oriented so that each column is a variable and each row a new measurement
#'
#' @return list containing the COVMTX, avg and tlen
#'
#' @seealso  \code{\link{lcovCreate}}
#' @export
#' @keywords internal
###################################################################################
lcovUpdate <- function (lcov, DATA){
	lcov$COVMTX= lcov$COVMTX +t(DATA)%*%DATA;
	lcov$avg=lcov$avg+colSums(DATA); #apply(DATA,2,sum);  #<--- better performance than before
	lcov$tlen=lcov$tlen+customSize(DATA,1);
	return(lcov)
}

###################################################################################
#' Fix a covariance object
#' 
#'   Computes the definitive covariance matrix
#'   and the average of the covariance object referenced by lcov
#'   after a series of update operations.
#'
#' @param lcov 			A list that contains all information about the handled covariance-structure
#'
#' @return returns the fixed covariance list structure
#'
#' @seealso  \code{\link{lcovCreate}} \code{\link{lcovPca}}
#' @export
#' @keywords internal
###################################################################################
lcovFix <- function (lcov){
	tlen=lcov$tlen;
	lcov$COVMTX= lcov$COVMTX/(tlen-1) - t(lcov$avg)%*%lcov$avg/tlen/(tlen-1);  
	lcov$avg=lcov$avg/tlen;
	return(lcov)
}

###################################################################################
#' Transform a covariance object
#' 
#'   Computes the definitive covariance matrix
#'   and the average of the covariance object referenced by lcov
#'   after a series of update operations.
#'
#' @param lcov 			A list that contains all information about the handled covariance-structure to be transformed
#' @param A				linear function by which covariance object is to be transformed
#'
#' @note lcovFix(lcov) has to be used before this function is applied
#'
#' @return returns the fixed covariance list structure
#'
#' @seealso  \code{\link{lcovFix}}
#' @export
#' @keywords internal
###################################################################################
lcovTransform <- function (lcov, A){ #TODO: untested
	lcov$avg<-lcov$avg %*% t(A)
	lcov$COVMTX<-lcov$COVMTX %*% t(A)  
	return(lcov)
}

###################################################################################
#' Principal Component Analysis on a covariance object
#' 
#'   Performs PCA _and_ whitening
#'   on the covariance object referenced by lcov.
#'   CAUTION: can be numerically instable if covariance matrix is singular,
#'   better use LCOV_PCA2 instead /W. Konen/
#'
#' @param lcov 			A list that contains all information about the handled covariance-structure
#' @param dimRange		A number or vector for dimensionality reduction: \cr 
#'         if it is a number: only the first components 1:dimRange are kept (those with largest eigenvalues)\cr            
#'         if it is a range: only the components in the range dimRange[1]..dimRange[2] are kept
#'
#' @note lcovFix(lcov) has to be used before this function is applied
#'
#' @return returns a list: 
#'   $W is the whitening matrix, $DW the dewhitening matrix and $D an array
#'   containing a list of the eigenvalues. $kvar contains the total
#'   variance kept in percent.
#'
#' @seealso  \code{\link{lcovFix}} \code{\link{lcovPca2}}
#' @export
#' @keywords internal
###################################################################################
lcovPca <- function (lcov,dimRange=NULL){
	if(is.null(dimRange)){
		dimRange=nrow(lcov$COVMTX)
	}
	#dimRange can be single value or two values that specify the interval
	dimInt=sfaGetIntRange(dimRange)
	
	svdResult=svd(lcov$COVMTX,nu=0,nv=ncol(lcov$COVMTX)) # singular value decomposition 
	
	D=svdResult$d;
	PC=svdResult$v;
	kvar = sum(D[dimInt])/sum(D);
	D=D[dimInt]^(-0.5); 
	#% CAUTION: this can lead to Inf, if any of 
    #% the diagonal elements of D is zero >>
    #% better use LCOV_PCA2 instead! /WK/
  
	Dmtx=diag(D); #construct a diagonal matrix from D values
	
	PC=t(PC[,dimInt]);
	W=Dmtx %*% PC; 
	DW=t(PC) %*% solve(Dmtx);#R solve() == inv() in matlab
	D=D^(-2);
	
	resList=list(D=D, W=W, DW=DW, kvar=kvar);
	return(resList)
}

###################################################################################
#'  Improved Principal Component Analysis on a covariance object
#' 
#'   Performs PCA _and_ whitening
#'   on the covariance object referenced by lcov. \cr \cr
#'   Difference to LCOV_PCA: null the rows of W (columns of DW) where the 
#'   corresponding eigenvalue in D is close to zero (more precisely: if
#'   lam/lam_max < EPS = 1e-7). This is numerically stable in the case where
#'   the covariance matrix is singular.\cr
#'  - Author: Wolfgang Konen, Cologne Univ., May'2009
#'
#' @param lcov 			A list that contains all information about the handled covariance-structure
#' @param dimRange		A number or vector for dimensionality reduction: \cr 
#'         if it is a number: only the first components 1:dimRange are kept (those with largest eigenvalues)\cr            
#'         if it is a range: only the components in the range dimRange[1]..dimRange[2] are kept
#'
#' @note lcovFix(lcov) has to be used before this function is applied
#'
#' @return returns a list: 
#'   $W is the whitening matrix, $DW the dewhitening matrix and $D an array
#'   containing a list of the eigenvalues. $kvar contains the total
#'   variance kept in percent.
#'
#' @seealso  \code{\link{lcovFix}} \code{\link{lcovPca}}
#' @export
#' @keywords internal
###################################################################################
lcovPca2 <- function (lcov,dimRange=NULL){
	if(is.null(dimRange)){
		dimRange=nrow(lcov$COVMTX)
	}
	#dimRange can be single value or two values that specify the interval
	dimInt=sfaGetIntRange(dimRange) 
	EPS=1e-7;
	svdResult=svd(lcov$COVMTX,nu=0,nv=ncol(lcov$COVMTX)) # singular value decomposition   
	
	d=svdResult$d;
	lammax=d[1]	
	PC=svdResult$v;
	kvar = sum(d[dimInt])/sum(d);#sum(d[dimInt],na.rm=TRUE)/sum(d,na.rm=TRUE); #TODO discuss: in one instance 1 of 299 elements of d was NaN
	
	d=d[dimInt];
	Dmtx=d^(-0.5);
	Dmtx[d<lammax*EPS]=0.0; 	#set very small values to zero to prevent numerical error amplification
	#Dmtx[is.na(d)]=0.0		#TODO discuss: also set values returned as NaN to zero=
	Dmtx=diag(Dmtx); #construct a diagonal matrix from D values

	PC=t(PC[,dimInt]);
	W=Dmtx %*% PC;
	DW=t(PC) %*% diag(d);
	
	resList=list(D=d, W=W, DW=DW, kvar=kvar);
	return(resList)
}