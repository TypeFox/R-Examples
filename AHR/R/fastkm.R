#' Fast Kaplan-Meier estimator
#'
#' This function calculates the Kaplan-Meier estimator for right-censored survival data,
#' at arbitrary time points. It can handle left-truncated and/or right-censored data with ties.
#' Avoids the overhead of the \code{survfit} or \code{prodlim} functions by stripping
#' away most of the features not needed here.
#'
#' @title fastkm
#' @param time vector of right-censored survival times
#' @param status censoring indicator for each element of \code{time} (0 = right-censored, 1 = event)
#' @param ltrunc vector of left-truncation times
#' @param left.limit indicates wether estimated survival function is left continuous
#' @param eval points at which the estimated survival function should be evaluated
#' @return A list containing the vectors \code{time}, \code{surv} and \code{variance}, and \code{n.atrisk}
#' @seealso \code{\link{survfit}} and \code{\link[prodlim]{prodlim}}
#' @export
#' @examples
#' T <- rexp(100)
#' C <- rexp(100)
#' Y <- pmin(T, C)
#' D <- T <= C
#' sort(fastkm(Y, D)$surv, decreasing=TRUE)
#' # should be exactly the same as
#' fit <- survfit(Surv(Y, D) ~ 1)
#' f <- approxfun(fit$time, fit$surv, f=0, rule=2, yleft=1)
#' f(fit$time)
fastkm <- function(time, status, ltrunc=rep.int(0, length(time)), left.limit=FALSE, eval=time) {    
    df <- data.frame(time=time, status=status, ltrunc=ltrunc)
    df <- df[order(df$time), ]
    .Call("R_fastkm", df$time, df$status, df$ltrunc, left.limit, eval, PACKAGE="AHR")
}
