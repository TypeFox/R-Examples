#' Summarize Rank-Based Linear Model Fits
#' 
#' Provides a summary similar to the traditional least squares fit.
#' 
#' 
#' @param object an object of class 'rfit', usually, a result of a call to
#' 'rfit'
#' @param \dots additional arguments
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' data(baseball)
#' fit<-rfit(weight~height,data=baseball)
#' summary(fit)
#' 
#' @export summary.rfit
summary.rfit <- function (object,...) {

  tauhat <- object$tauhat
  n<-length(object$y)
  pp1 <- object$qrx1$rank
  est <- object$coef
  ses <- sqrt(diag(vcov(object)))
  tstat <- est/ses
  pval <- 2 * pt(-abs(tstat), n - pp1)
  coef <- cbind(est, ses, tstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "t.value","p.value")
  dt <- drop.test(object)
  R2 <- (dt$df1/dt$df2 * dt$F)/(1 + dt$df1/dt$df2 * dt$F)
  ans <- list(coefficients = coef, dropstat = dt$F, droppval = dt$p.value, 
    R2 = R2)
  ans$call <- object$call
  class(ans) <- "summary.rfit"
  ans

}
