# This code is based on the matlab packages SFA-tk 1.0 and the expansion V2.X

# Package Description for Roxygene:
#' Slow Feature Analysis in R
#'
#' \tabular{ll}{
#' Package: \tab rSFA\cr
#' Type: \tab Package\cr
#' Version: \tab 1.04\cr
#' Date: \tab 17.12.2014\cr
#' Maintainer: \tab Martin Zaefferer \email{martin.zaefferer@@gmx.de}\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Slow Feature Analysis in R, ported to R based on the matlab versions SFA toolkit 1.0 by Pietro Berkes and SFA
#'    toolkit 2.8 by Wolfgang Konen.
#'
#' @title Slow Feature Analysis in R
#' @name rSFA-package
#' @aliases rSFA
#' @docType package
#' @author Wolfgang Konen \email{wolfgang.konen@@fh-koeln.de}, Martin Zaefferer, Patrick Koch; 
#'			Bug hunting and testing by Ayodele Fasika, Ashwin Kumar, Prawyn Jebakumar
#' @keywords slow feature analysis timeseries classification
#' @references \url{http://gociop.de/research-projects/sfa/}
# @useDynLib rSFA  #remark: deprecated! not working properly, thus removed
#End of Package Description
NA #NULL, ends description without hiding first function

###################################################################################
#' The SFA1 algorithm, linear SFA. 
#' 
#'    Y = sfa1(X) performs linear Slow Feature Analysis on the input data
#'    X and returns the output signals Y ordered by increasing temporal
#'    variation, i.e. the first signal Y[,1] is the slowest varying one,
#'    Y[,2] the next slowest and so on. The input data have to be organized 
#'    with each variable in a column and each data (time) point in a
#'    row, i.e. X(t,i) is the value of variable nr. i at time t.
#'
#' @param x 			Input data, each column a different variable
#'
#' @return list \code{sfaList} with all learned information, where \code{sfaList$y} contains the outputs
#'
#' @seealso  \code{\link{sfaStep}} \code{\link{sfa1Create}} \code{\link{sfaExecute}}
#' @export
###################################################################################
sfa1 <- function (x){
	#number of input signals
	n=ncol(x)
	
	#create a SFA list
	sfaList=sfa1Create(n)

	#perform the preprocessing step
	sfaList=sfaStep(sfaList, x, "preprocessing");

	#close the algorithm
	sfaList=sfaStep(sfaList, NULL, 'sfa');

	#compute the output signal
	y = sfaExecute(sfaList, x);

	sfaList$y=y
	return(sfaList)
}

###################################################################################
#' The SFA2 algorithm, SFA with degree 2 expansion. 
#' 
#'    Y = sfa2(X) performs expanded Slow Feature Analysis on the input data
#'    X and returns the output signals Y ordered by increasing temporal
#'    variation, i.e. the first signal Y[,1] is the slowest varying one,
#'    Y[,2] the next slowest varying one and so on. The input data have to
#'    be organized with each variable in a column and each data (time) point in a
#'    row, i.e. X(t,i) is the value of variable i at time t.
#'	  By default an expansion to the space of 2nd degree polynomials is done,
#'	  this can be changed by using different functions for xpDimFun and sfaExpandFun.
#'
#' @param x 				input data
#' @param method 		eigenvector calculation method:	="SVDSFA" for singular value decomposition (recommended) or 
#' 						 		  ="GENEIG" for generalized eigenvalues (unstable!).	GENEIG is not implemented in the current version, since
#' 						  		R lacks an easy option to calculate generalized eigenvalues.
#' @param ppType		preprocessing type: ="PCA" (principal component analysis) or ="SFA1" (linear sfa)
#' @param xpDimFun			 function to calculate dimension of expanded data
#' @param sfaExpandFun	 function to expand data 
#'
#' @return list \code{sfaList} with all SFA information, among them are
#'    \item{\code{y}}{ a matrix containing  the output Y (as described above) }
#'    \item{-}{ all input parameters to \code{\link{sfa2Create}}  }
#'    \item{-}{ all elements of \code{sfaList}  as specified in \code{\link{sfa2Step}}}
#'
#' @seealso  \code{\link{sfa2Step}} \code{\link{sfa2Create}} \code{\link{sfaExecute}} \code{\link{sfa1}} 
#' @examples
#' ## prepare input data for simple demo
#' t=seq.int(from=0,by=0.011,to=2*pi)
#' x1=sin(t)+cos(11*t)^2
#' x2=cos(11*t)
#' x=data.frame(x1,x2)
#' ## perform sfa2 algorithm with data
#' res = sfa2(x)
#' ## plot slowest varying function of result
#' plot(t, res$y[,1],type="l",main="output of the slowest varying function")
#' ## see http://www.scholarpedia.org/article/Slow_feature_analysis#The_algorithm
#' ## for detailed description of this example
#' @export
###################################################################################
sfa2 <- function (x,method="SVDSFA",ppType="PCA", xpDimFun=xpDim, sfaExpandFun=sfaExpand){
	#number of input signals
	n=ncol(x)
		
	#create a SFA list
	sfaList=sfa2Create(n, xpDimFun(n), ppType, xpDimFun=xpDimFun, sfaExpandFun=sfaExpandFun)

	#perform the preprocessing step
	sfaList=sfaStep(sfaList, x, "preprocessing");

	#perform expansion step
	sfaList=sfaStep(sfaList, x, "expansion");

	#close the algorithm
	sfaList=sfaStep(sfaList, NULL, "sfa", method); #WHY WAS HERE second argument originally NULL and now it is x??? 1.0 -> 2.x matlab version

	#compute the output signal
	y = sfaExecute(sfaList, x);

	# check unit variance constraint of slow (training) signals
	# /WK/08/2009/
	Cslow = diag(t(y)%*%y)/customSize(y,1); #TODO check if it is necessary (maybe only for debugging? performance issue?)
	if (max(abs(Cslow-1))>0.1 & method!="SVDSFA"){
		warning("
			It is recommended to use method=SVDSFA in SFA-step.
			Some of the slow signals y_i have a variance <> 1:",
			paste(Cslow," "));
	}	
	sfaList$y=y
	return(sfaList)
}