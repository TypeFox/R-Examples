#' Summarize a fitted \code{grpreg} object
#' @description Summarizes a fitted penalized regression model with '\code{grpreg}' class.
#'
#' @param object A "\code{summary.grpreg}" object.
#' @param lambda A regularization parameter at which to summarize.
#' @param digits Number of digits past the decimal point to print out. The default is \code{4}.
#' @param ... Optional arguments passed to other methods.
#'
#' @details The default \code{lambda} is the one at which we obtain the minimum loss
#' (i.e., negative log-likelihood value), if \code{lambda} is not supplied.
#' @return A list with class "\code{summary.grpreg}" containing the following components:
#' \item{n}{Number of observations.}
#' \item{p}{Number of screened predictors.}
#' \item{penalty}{The penalty applied to the model.}
#' \item{model}{The type of model.}
#' \item{family}{The link function.}
#' \item{criterion}{The screening criterion.}
#' \item{lambda}{The default or specified regularization parameter.}
#' \item{beta}{The estimates of coefficients at the specified \code{lambda}.}
#' \item{iter}{The number of iterations at the specified \code{lambda}.}
#' \item{df}{The estimates of effective number of model parameters at the specified
#'  \code{lambda}.}
#' \item{call}{The function call.}
#' Additional elements are contained for the case in which \code{family = "gaussian"}:
#' \item{r.squared}{The r.squared.}
#' \item{snr}{The signal-to-noise ratio.}
#' \item{scale}{The scale parameter estimate.}
#' and following elements for the case in which \code{family = "poisson"} or \code{"binomial"}:
#' \item{logLik}{The negative log-likelihood values for the fitted model.}
#' \item{aic}{Akaike's information criterion (AIC).}
#' \item{bic}{Bayesian information criterion (BIC).}
#' \item{aicc}{The AIC with a correction for finite sample sizes (AICC).}
#' \item{pe}{The prediction error for \code{family = "binomial"}.}
#'
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}, \code{\link{summary.cv.grpreg}}
#'
#' @importFrom stats var
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats logLik
#'
#' @export
summary.grpreg <- function(object, lambda = NULL, digits = 4, ...) {
  Min <- ifelse(is.null(lambda),which.min(object$loss),
                which.min(abs(object$lambda - lambda)))
  if (Min > length(object$lambda))
    stop(paste("please specify 'lambda' between",min(object$lambda),
               "and",max(object$lambda)))
  beta <- object$beta[,Min]
  lambda <- ifelse(is.null(lambda),object$lambda[Min],lambda)
  callArg <- object$call
  if (is.null(callArg$formula))
    y <- eval(callArg$y)
  else {
    data <- eval(callArg$data)
    yresp <- as.character(callArg$formula)[2]
    y <- data[,yresp]
  }
  p <- length(beta) - 1
  n <- length(y)
  LL <- logLik(object)
  ypred <- predict(object)
  if (object$family == "gaussian") {
    R.squared <- 1 - sum((y - ypred)^2)/sum((y - mean(y))^2)
    snr <- var(mean(y) - ypred)/var(y - ypred)
    result <- list(beta = beta, r.squared = R.squared, snr = snr, lambda = lambda,
                   scale = sd(y - ypred))
  }
  else {
    aic <- AIC(object)[Min]
    bic <- BIC(object)[Min]
    aicc <- aic + 2*p*(p - 1)/(n - p - 1)
    result <- list(beta = beta, logLik = LL[Min], aic = aic, bic = bic, aicc = aicc,
                   lambda = lambda)
    if (object$family == "binomial") {
      tab <- table(y,ypred)
      result <- c(result, list(pe = 1 - sum(diag(tab))/sum(tab)))
    }
  }
  penal <- switch(object$penalty,grSCAD = "Group SCAD",
                  grLasso = "Group Lasso", grMCP = "Group MCP",
                  gel = "Group exponential",cMCP = "Composite MCP")
  criteria <- switch(object$criterion, gSIS = "SIS", gHOLP = "HOLP",
                     gAR2 = "AR2", gDC = "DC")
  model <- switch(object$family, gaussian = "linear",binomial = "logistic",
                  poisson = "poisson")
  result <- c(list(n = n, p = p, penalty = penal, model = model,
                   family = object$family, criteria = criteria), result,
              list(iter = object$iter[Min], df = object$df[Min],call = object$call))
  class(result) <- "summary.grpreg"
  return(result)
}
