#' Compute standardized beta coeffizients for linear regression models
#' @export
#' @param MOD A model of class \code{lm}.
#' @return A data.frame with two columns
#' \tabular{ll}{
#'   \code{beta} \tab standardized beta coefficients\cr
#'   \code{se.beta} \tab standard errors for the beta coefficients\cr
#' }
#' @examples
#' data <- pisa2012che
#'
#' ## linear regression models with numerical covariates only
#' fit1 <- lm(MATH ~ ESCS + USEMATH, data)
#' betas.lm(fit1)
#'
#' ## ...and with interaction terms
#' fit1.1 <- lm(MATH ~ ESCS * USEMATH, data)
#' betas.lm(fit1.1)
#'
#' ## linear regression models with numerical and factorial covariates
#' fit2 <- lm(MATH ~ ESCS + USEMATH + ST04Q01 + FAMSTRUC + ST28Q01, data)
#' betas.lm(fit2)
#'
#' ## ...and with interaction terms
#' fit2.1 <- lm(MATH ~ ESCS + USEMATH + ST04Q01 + FAMSTRUC * ST28Q01, data)
#' betas.lm(fit2.1)
#'
#' ## weighted linear regression models
#' fit3 <- lm(MATH ~ ESCS + USEMATH, data, weights = W_FSTUWT)
#' betas.lm(fit3)
#'
#' fit4 <- lm(MATH ~ ESCS + USEMATH + ST04Q01 + FAMSTRUC + ST28Q01, data, weights = W_FSTUWT)
#' betas.lm(fit4)
#'
#' ## ...with interaction terms
#' fit3.1 <- lm(MATH ~ ESCS * USEMATH, data, weights = W_FSTUWT)
#' betas.lm(fit3.1)
#'
#' fit4.1 <- lm(MATH ~ ESCS + USEMATH + ST04Q01 + FAMSTRUC * ST28Q01, data, weights = W_FSTUWT)
#' betas.lm(fit4.1)

betas.lm <- function (MOD) {
  if(class(MOD) != "lm")
    stop("Object must be of class 'lm'.")

  if(is.null(MOD$qr))
    stop("Please refit model with 'qr = TRUE'.")

  ## coefficients w/o intercept
  b <- MOD$coefficients[-1]

  ## stand. errors w/o intercept
  se <- summary(MOD)$coefficients[-1,2]

  ## compute sd w/ and w/o weights
  X <- qr.X(MOD$qr)
  if(is.null(w <- MOD$weights)) {
    sdx <- apply(X, 2, sd)[-1]
  } else {
    sdx <- apply(X, 2, function(Y) sqrt(my.wtd.var(Y, w)))[-1]
  }

  ## sd of response variable
  sdy <- sd(MOD$model[,1])

  ## beta = b * sd(x)/sd(y)
  beta <- b * sdx/sdy
  se.b <- se * sdx/sdy

  return(data.frame(beta=beta, se.beta=se.b))
}

#' Compute standardized beta coeffizients for robust linear regression models
#' @export
#' @param object A model of class \code{lmRob}.
#' @param classic Logical TRUE for classic covariance estimation.
#' @return Vector with standardized beta coefficients.
#' @importFrom robust covRob
#' @importFrom robust covClassic
#' @examples
#' library(robust)
#' data <- pisa2012che
#'
#' ## robust estimation of betas
#' fit1 <- lmRob(MATH ~ ESCS, data)
#' betas.lmr(fit1)
#'
#' ## example where robust variance cannot be computed,
#' ## instead the classical variance is used.
#' fit2 <- lmRob(MATH ~ ESCS + USEMATH, data)
#' betas.lmr(fit2)

betas.lmr <- function (object, classic = FALSE) {
  if(class(object) != "lmRob")
    stop("Object must be of class 'lmRob'")
  model <- object$model
  #num <- sapply(model, is.numeric)  # numeric vars only
  b <- object$coefficients[-1]  # final coefficients w/o intercept
  ## compute robust covariance
  covr <- NULL
  try(covr <- diag(covRob(model)$cov), silent = TRUE)
  if(is.null(covr) & classic == FALSE)
    warning("covRob() coud not be computed, instead covClassic() was applied.")
  ## compute classic covariance if robust failed
  if(is.null(covr) | classic == TRUE)
    covr <- diag(covClassic(sapply(model, as.numeric))$cov)
  sx <- sqrt(covr[-1])  # standard deviation of x's
  sy <- sqrt(covr[1])  # standard deviation of y
  beta <- b * sx/sy
  return(beta)
}

# form package Hmisc
# @param x a numeric vector
# @param weights a numeric vector of weights
# @param normwt specify normwt=TRUE to make weights sum to length(x) after deletion of NAs
# @param na.rm  set to FALSE to suppress checking for NAs
# @return weighted variance.

my.wtd.var <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE)
{
  if (!length(weights)) {
    if (na.rm)
      x <- x[!is.na(x)]
    return(var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt)
    weights <- weights * length(x)/sum(weights)
  sw <- sum(weights)
  xbar <- sum(weights * x)/sw
  sum(weights * ((x - xbar)^2))/(sw - sum(weights^2)/sw)
}
