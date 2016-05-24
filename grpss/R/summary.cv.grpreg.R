#' Summarize a fitted \code{cv.grpreg} object
#' @description Summarizes a fitted penalized regression model with '\code{cv.grpreg}' class.
#'
#' @param object A fitted "\code{cv.grpreg}" object.
#' @param digits Number of digits past the decimal point to print out. The default is \code{4}.
#' @param ... Optional arguments passed to other methods.
#' @details This function is similar to \code{summary.cv.grpreg} in \code{grpreg} package, but
#' gives several different results. The \code{lambda} here is only valid for "\code{grpreg}"
#' object. See \code{\link{summary.grpreg}}.
#' @return A list with class "\code{summary.cv.grpreg}" containing the following components:
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
#' \code{lambda}.}
#' \item{cve}{The average cross-validation error at specified \code{lambda}.}
#' \item{call}{The function call.}
#' Additional elements are contained for the case in which \code{family = "gaussian"}:
#' \item{r.squared}{The r.squared.}
#' \item{snr}{The signal-to-noise ratio.}
#' \item{scale}{The scale parameter estimate (sigma).}
#' and following elements for the case in which \code{family = "poisson"} or \code{"binomial"}:
#' \item{logLik}{The negative log-likelihood values for the fitted model.}
#' \item{aic}{Akaike's information criterion (AIC).}
#' \item{bic}{Bayesian information criterion (BIC).}
#' \item{aicc}{The AIC with a correction for finite sample sizes (AICC).}
#' \item{pe}{The prediction error for \code{family = "binomial"}.}
#'
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}, \code{\link{summary.grpreg}}
#' @export
#' @importFrom stats predict

summary.cv.grpreg <- function(object,digits = 4,...) {
  callArg <- object$call
  if (is.null(callArg$formula))
    y <- eval(callArg$y)
  else {
    data <- eval(callArg$data)
    yresp <- as.character(callArg$formula)[2]
    y <- data[,yresp]
  }
  ypred <- predict(object)
  index <- which.min(abs(object$lambda - object$lambda.min))
  beta <- object$fit$beta[,index]
  p  <- length(beta) - 1
  n <- length(y)
  LL <- object$fit$loss[index]
  if (object$fit$family == "gaussian") {
    R.squared <- 1 - sum((y - ypred)^2)/sum((y - mean(y))^2)
    snr <- var(mean(y) - ypred)/var(y - ypred)
    result <- list(beta = beta, r.squared = R.squared, snr = snr,
                   scale = sd(y - ypred), lambda = object$lambda.min)
  }
  else {
    index <- which.min(abs(object$lambda - object$lambda.min))
    aic <- 2*LL + 2*p
    aicc <- aic + 2*p*(p + 1)/(n - p - 1)
    bic <- 2*LL + p*log(n)
    result <- list(beta = beta, logLik = LL, aic = aic, bic = bic, aicc = aicc,
                   lambda = object$lambda.min)
    if (object$fit$family == "binomial") {
      tab <- table(y,ypred)
      result <- c(result, list(pe = 1 - sum(diag(tab))/sum(tab)))
    }
  }
  penal <- switch(object$fit$penalty,grSCAD = "Group SCAD",
                  grLasso = "Group Lasso", grMCP = "Group MCP",
                  gel = "Group exponential",cMCP = "Composite MCP")
  criteria <- switch(object$criterion, gSIS = "SIS", gHOLP = "HOLP",
                     gAR2 = "AR2", gDC = "DC")
  model <- switch(object$fit$family, gaussian = "linear",
                  binomial = "logistic",poisson = "poisson")
  result <- c(list(n = n, p = p, penalty = penal, model = model,
                   family = object$fit$family, criteria = criteria), result,
              list(cve = object$cve[index],iter = object$iter[index], df = object$df[index],
                   call = object$call))
  class(result) <- "summary.cv.grpreg"
  return(result)
}
