# TODO: Add comment
# 
# Author: soueidat
###############################################################################

###################################################################################
##' Constructor of Input class
##'
##' This class contains the input paramaters need to run the algorithm.
##'
##' \describe{
##'   \item{coefs}{In the univariate case this will be the transpose of the coefficients matrix. In the multivarite case this matrix will be the concatenation of the coefs matrix for each dimension of the multivariate functinal object.}
##'   \item{basisProd}{In the univariate case this is the matrix of the inner product between the basis functions.In the multivariate case the m_basisProd member will be a block diagonal matrix and each block will be the matrix of the inner product between the basis functions of each dimension of the data.} 
##'   \item{K}{the number of clusters.}
##'   \item{thd}{the threshold in the scree test to select the dimensions of the curves.}
##'   \item{increaseDimension}{A logicla paramater, if true the dimensions will be constraint to increase after each iteration. A false mean that the dimensions will take their values according to the scree test.}
##'   \item{hard}{A logical parameter, if true we initialize randomly the model with "hard" weights  A logical parameter, if true we initialize randomly the model with hard weights  (weights of each curves taking 0 or 1 according to the class membership of the curves). if false we initialize randomly with "soft" weights (weights of each curves taking a probabilities according to the class membership of the curves).}
##'   \item{fixedDimension}{A vector of size "K" which contains the dimensions in the case of running algorithm with fixed dimensions.}  
##'   \item{epsilon}{The stoping criterion, we stop run the algorithm if the difference between two successive loglikelihood is less than epsilon.}
##'   \item{nbInit}{The number of initialization to be achieve, befor running the long algorithm.}
##'   \item{nbIterInit}{The maximum number of iterations in each initialization.}
##'   \item{nbIteration}{The maximum number of iteration in the long algorithm.}
##'
##' }
##'
##'
##' @name Input-class
##' @rdname Input-class
## @exportClass Input
##'

setClass(
		Class = "Input",
		representation = representation(
				coefs = "matrix",
				basisProd = "matrix",
				K = "numeric",
				thd= "numeric",
				increaseDimension="logical",
				hard="logical",
				fixedDimension="integer",
				epsilon = "numeric",
				nbInit = "numeric",
				nbIterInit= "numeric",
				nbIteration= "numeric"
		),
		prototype = prototype(
				coefs = matrix(nrow=0,ncol=0),
				basisProd = matrix(nrow=0,ncol=0),
				K = numeric(0),
				thd= numeric(0),
				increaseDimension=FALSE,
				hard=TRUE,
				fixedDimension=integer(0),
				epsilon = numeric(0),
				nbInit = numeric(0),
				nbIterInit= numeric(0),
				nbIteration= numeric(0)
		)
)
