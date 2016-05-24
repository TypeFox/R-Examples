##' Constructor of Output class
##'
##' This class contains the parameters in the output after running classification.
##'
##' \describe{
##'   \item{tik}{a matrix of size (number of Curves) x (K), each column contains the weights of the curves in the corresponding class.}
##'   \item{cls}{a vector of size number of curves, containing the index of the class for each curve.}
##'   \item{proportion}{a matrix of size 1xnbClust (number of clusters), containing the estimated mixture proportions.}
##'   \item{loglikelihood}{the estimated log-likelihood.}
##'   \item{aic}{the value of AIC criterion.}
##'   \item{bic}{the value of BIC criterion.}
##'   \item{icl}{the value of ICL criterion.}
##'   \item{dimensions}{a vector of size nbClust of the dimensions of the specifique dimensions of the functional data in each class.}
##'   \item{dimTotal}{a matrix of size nbClust x nbRunIteration, where nbRunIteration is the number of iterations before the algorithm converge. Each column of the dimTotal matrix contain the dimensions on the coresponding iteration.}
##'   \item{V}{principal components variances per cluster}
##'   \item{empty}{logical parameter, and empty=TRUE if we have an empty class}
##'
##' }
##'
##'
##' @name Output-class
##' @rdname Output-class
## @exportClass Output
##'

setClass(
		Class="Output",
		representation = representation(
				tik="matrix",
				cls= "numeric",
				proportions = "matrix",
				loglikelihood = "numeric",
				loglikTotal="numeric",
				aic="numeric",
				bic="numeric",
				icl="numeric",
				dimensions="numeric",
				dimTotal="matrix",
				V="matrix",
				empty="logical"
		),
		prototype = prototype(
				tik=matrix(nrow=0,ncol=0),
				cls= numeric(0),
				proportions = matrix(nrow=0,ncol=0),
				loglikelihood = numeric(0),
				loglikTotal=numeric(0),
				aic=numeric(0),
				bic=numeric(0),
				icl=numeric(0),
				dimensions=numeric(0),
				dimTotal=matrix(nrow=0,ncol=0),
				V=matrix(nrow=0,ncol=0),
				empty=FALSE
		)
)
