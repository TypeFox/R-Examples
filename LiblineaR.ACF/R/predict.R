### Documentation ####
#' Predictions with LiblineaR.ACF model
#' 
#' The function applies a classification model produced by the \code{LiblineaR.ACF} function to every row of a
#' data matrix and returns the model predictions.
#' 
#' @param object Object of class \code{"LiblineaR.ACF"}, created by
#'   \code{LiblineaR.ACF}.
#' @param newx An n x p matrix containing the new input data. A vector will be
#'   transformed to a n x 1 matrix. A sparse matrix (from SparseM package) will
#'   also work.
#' @param decisionValues Logical indicating whether model decision values should
#'   be computed and returned. Default is \code{FALSE}.
#' @param ... Currently not used
#' 
#' @return 	By default, the returned value is a list with a single entry:
#' \item{predictions}{A vector of predicted labels.}
#' If \code{decisionValues} is set to \code{TRUE}, an additional entry is returned:
#' \item{decisionValues}{An n x k matrix (k number of classes) of the model
#'   decision values. The columns of this matrix are named after class labels.}
#' 
#' @references
#' 	\itemize{
#' \item 
#' For more information on LIBLINEAR itself, refer to:\cr
#' R.-E. Fan, K.-W. Chang, C.-J. Hsieh, X.-R. Wang, and C.-J. Lin.\cr
#' \emph{LIBLINEAR: A Library for Large Linear Classification,}\cr
#' Journal of Machine Learning Research 9(2008), 1871-1874.\cr
#' \url{http://www.csie.ntu.edu.tw/~cjlin/liblinear}
#' }
#' 
#' @author Thibault Helleputte \email{thibault.helleputte@@dnalytics.com} and 
#'   Pierre Gramme \email{pierre.gramme@@dnalytics.com}.\cr 
#'   Modified by Aydin Demircioglu.\cr
#'   Based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin
#' 
#' @note If the data on which the model has been fitted have been centered
#'   and/or scaled, it is very important to apply the same process on the
#'   \code{newx} data as well, with the scale and center values of the training
#'   data.
#' 
#' @seealso \code{\link{LiblineaR.ACF}}
#' 
#' @keywords classif multivariate models optimize classes
#'
#' @export

### Implementation ####
predict.LiblineaR.ACF<-function(object, newx, decisionValues=FALSE,...){

	# <Arg preparation>
	
	error=c()
	
  if(sparse <- inherits(newx, "matrix.csr")){
    if(requireNamespace("SparseM",quietly=TRUE)){
		  # trying to handle the sparse martix case
		  newx = SparseM::t(SparseM::t(newx)) # make sure column index are sorted
		  n = newx@dimension[1]
		  p = newx@dimension[2]
    } else {
      stop("newx inherits from 'matrix.csr', but 'SparseM' package is not available. Cannot proceed further. Use non-sparse matrix or install SparseM.")
    }
	} else {
		# Nb samples
		n=dim(newx)[1]
		# Nb features
		p=dim(newx)[2]
	}
	
	# Type 
	if(!object$Type %in% c(0:7,11:13)){
		stop("Invalid model object: Wrong value for 'type'. Must be an integer between 0 and 7  or between 11 and 13 included.\n")
	}
	
	# Bias
	b = if(object$Bias) 1 else -1
	
	# Returned probabilities default storage preparation
	Probabilities=matrix(data=-1)

	# Returned labels storage preparation
	Y=matrix(ncol=n,nrow=1,data=0)
	
	# Returned decision values default storage preparation
	DecisionValues=matrix(data=-1)

	# Returned decision values storage preparation
	if(decisionValues) {
        DecisionValues=matrix(ncol=n*length(object$ClassNames),nrow=1,data=0)
	}
	
	# Codebook for labels
	cn=c(1:length(object$ClassNames))
	
	#
	# </Arg preparation>
	
	# as.double(t(X)) corresponds to rewrite X as a nxp-long vector instead of a n-rows and p-cols matrix. Rows of X are appended one at a time.
	# as we only deal with classification, we put probA always to false
	proba = FALSE
	
	ret <- .C(
		"predictLinear",
		as.double(Y),
		as.double(if(sparse) newx@ra else t(newx)),
		as.double(object$W),
		as.integer(decisionValues),
		as.double(t(DecisionValues)),
		as.integer(proba),
		as.double(t(Probabilities)),
		as.integer(object$NbClass),
		as.integer(p),
		as.integer(n),
		# sparse index info
		as.integer(sparse),
		as.integer(if(sparse) newx@ia else 0),
		as.integer(if(sparse) newx@ja else 0),

		as.double(b),
		as.integer(cn),
		as.integer(object$Type),
		PACKAGE="LiblineaR.ACF"
		)
	
	result=list()
    result$predictions=object$ClassNames[as.integer(ret[[1]])]
	
	if(decisionValues){
		result$decisionValues=matrix(ncol=length(object$ClassNames),nrow=n,data=ret[[5]],byrow=TRUE)
		colnames(result$decisionValues)=object$ClassNames
	}

	return(result)

}
