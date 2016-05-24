#' Studentized Residuals for Rank-Based Regression
#' 
#' Returns the Studentized residuals based on rank-based estimation.
#' 
#' 
#' @param model an object of class rfit
#' @param \dots additional arguments. currently not used.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' x<-runif(47)
#' y<-rcauchy(47)
#' qqnorm(rstudent(fit<-rfit(y~x)))
#' plot(x,rstudent(fit)) ; abline(h=c(-2,2))
#' 
#' @export rstudent.rfit
"rstudent.rfit" <- function (model,...) {
  fit<-model
  ehat <- fit$resid
  n <- length(ehat)
  p <- fit$qrx1$rank-1
  sigmahat <- mad(ehat)
  deltas <- sum(abs(ehat))/(n - p)
  delta <- disp(fit$betahat, fit$x, fit$y, fit$scores)/(n - p)
  k2 <- (fit$tauhat/sigmahat)^2 * (2 * delta/fit$tauhat - 1)
  if (fit$symmetric) {
    h <- hat(fit$x)
    s <- suppressWarnings(sigmahat * sqrt(1 - k2 * h))
    s[is.na(s)] <- sigmahat * sqrt(1 - h)[is.na(s)]
  } else {
    hc <- hat(as.matrix(qr.Q(fit$qrx1)[,2:(p+1)]), intercept=FALSE)
    k1 <- (fit$taushat/sigmahat)^2 * (2 * deltas/fit$taushat - 1)
    s <- suppressWarnings(sigmahat * sqrt(1 - k1/n - k2 * hc))
    s[is.na(s)] <- sigmahat * sqrt(1 - hc)[is.na(s)]
  }
  ehat/s
}
