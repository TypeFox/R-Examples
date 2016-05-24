#' Drop (Reduction) in Dispersion Test
#' 
#' Given two full model fits, this function performs a reduction in disperion
#' test.
#' 
#' Rank-based inference proceedure analogous to the traditional (LS) reduced
#' model test.
#' 
#' @param fitF An object of class rfit.  The full model fit.
#' @param fitR An object of class rfit.  The reduced model fit.
#' @return % ~Describe the value returned % If it is a LIST, use \item{F}{Value
#' of the F test statistic} \item{p.value}{The observed significance level of
#' the test (using an F quantile)} \item{RD}{Reduced model dispersion minus
#' Full model dispersion} \item{tauhat}{Estimate of the scale parameter (using
#' the full model residuals)} \item{df1}{numerator degrees of freedom}
#' \item{df2}{denominator degrees of freedom} %\item{comp1 }{Description of
#' 'comp1'} %\item{comp2 }{Description of 'comp2'} % ...
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' y<-rnorm(47)
#' x1<-rnorm(47)
#' x2<-rnorm(47)
#' fitF<-rfit(y~x1+x2)
#' fitR<-rfit(y~x1)
#' drop.test(fitF,fitR)
#' 
#' @export drop.test
drop.test <- function (fitF, fitR = NULL) {

  pp1 <- fitF$qrx1$rank

  if (is.null(fitR)) {
    rd <- disp(rep(0, ncol(fitF$x)), fitF$x, fitF$y, fitF$scores) - 
      disp(fitF$betahat, fitF$x, fitF$y, fitF$scores)
    df1 <- pp1 - 1
  } else {
    if( !all(abs( qr.fitted(fitF$qrx1,qr.Q(fitR$qrx1)) - qr.Q(fitR$qrx1) ) < .Machine$double.eps ^ 0.5 ) ) stop('Reduced model must be a subset of full model')
    rd <- disp(fitR$betahat, fitR$x, fitR$y, fitR$scores) - 
      disp(fitF$betahat, fitF$x, fitF$y, fitF$scores)
    df1 <- length(fitF$betahat) - length(fitR$betahat)
  }

  if( rd < 0 ) stop( "drop.test: negative reduction in dispersion found\n",
	"try starting full model at reduced model\n",
	"see help(drop.test) for more information" )

  df2 <- length(fitF$y) - pp1
  test <- (rd/df1)/(fitF$tauhat/2)
  pval <- 1 - pf(test, df1, df2)
  ans <- list(F = test, p.value = pval, RD = rd, tauhat = fitF$tauhat, 
    df1 = df1, df2 = df2)
  class(ans) <- "drop.test"
  ans

}
