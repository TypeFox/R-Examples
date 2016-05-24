#' Variance-Covariance Matrix for Rank-Based Regression
#' 
#' Returns the variance-covariance matrix of the regression estimates from an
#' object of type rfit.
#' 
#' 
#' @param object an object of type rfit
#' @param intercept logical. If TRUE include the variance-covariance estimates
#' corresponding to the intercept
#' @param \dots additional arguments
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @export vcov.rfit
vcov.rfit <- function (object, intercept = NULL, ...) {

  Q<-qr.Q(object$qrx1)
  q1<-Q[,1]
  Q2<-Q[,2:object$qrx1$rank]

  xxpxi<-object$x%*%chol2inv(chol(crossprod(object$x)))
  A1<-crossprod(q1,xxpxi) ; A2<-crossprod(Q2,xxpxi)
  object$taushat^2*crossprod(A1)+object$tauhat^2*crossprod(A2)

}
