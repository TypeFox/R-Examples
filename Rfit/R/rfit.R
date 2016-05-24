#' Rank-based Estimates of Regression Coefficients
#' 
#' Minimizes Jaeckel's dispersion function to obtain a rank-based solution for
#' linear models.
#' 
#' Rank-based estimation involves replacing the L2 norm of least squares
#' estimation with a pseudo-norm which is a function of the ranks of the
#' residuals. That is, in rank estimation, the usual notion of Euclidean
#' distance is replaced with another measure of distance which is referred to
#' as Jaeckel's (1972) dispersion function. Jaeckel's dispersion function
#' depends on a score function and a library of commonly used score functions
#' is included.  e.g. Wilcoxon and sign score functions. If an inital fit is
#' not supplied (i.e. yhat0 = NULL) then inital fit is based on an L1 fit via
#' rq.
#' 
#' @aliases rfit rfit.default
#' @param formula an object of class formula
#' @param data an optional data frame
#' @param subset an optional argument specifying the subset of observations to
#' be used
#' @param yhat0 an n by vector of initial fitted values, default is NULL
#' @param scores an object of class 'scores'
#' @param symmetric logical.  If 'FALSE' uses median of residuals as estimate
#' of intercept
#' @param TAU version of estimation routine for scale parameter.  F0 for
#' Fortran, R for (slower) R, N for none
#' @param \dots additional arguments to be passed to fitting routines
#' @return \item{coefficients}{estimated regression coefficents with intercept}
#' \item{residuals}{the residuals, i.e. y-yhat} \item{fitted.values}{ yhat = x
#' betahat} %\item{scores}{ score function used in estimation}
#' %\item{x}{original design matrix} %\item{y}{original response vector}
#' \item{xc}{centered design matrix} \item{tauhat}{estimated value of the scale
#' parameter tau} \item{taushat}{estimated value of the scale parameter tau_s}
#' %\item{symmetric}{ \item{betahat}{estimated regression coefficents}
#' \item{call}{Call to the function}
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso %~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{summary.rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449
#' - 1458.
#' 
#' Jureckova, J. (1971). Nonparametric estimate of regression coefficients.
#' \emph{Annals of Mathematical Statistics}, 42, 1328 - 1338.
#' @keywords nonparametric robust regression
#' @examples
#' 
#' data(baseball)
#' data(wscores)
#' fit<-rfit(weight~height,data=baseball)
#' summary(fit)
#' 
#' @export rfit
rfit <- function (formula, data=list(),...) 
UseMethod("rfit")
