#' Add All Possible Single Terms to a GWRM Model
#'
#' Compute all the single terms in the scope argument that can be added to the GWRM model, fit those models and compute a table of the changes in fit.
#'
#' @param object	a fitted object of class inheriting from \code{"gw"}.
#' @param scope a formula giving the terms to be considered for adding.
#' @param test \code{"none"}, which considers the AIC criterion, or \code{Chisq}, which is the likelihood-ratio test.
#' @param k the penalty constant in AIC / Cp.
#' @param trace	if \code{TRUE}, print out progress reports.
#' @param ...	further arguments passed to or from other methods.
#'
#' @return An object of class \code{"anova"} summarizing the differences in fit between the models.
#'
#' @importFrom stats pchisq add.scope update.formula extractAIC nobs formula update as.formula
#' @examples
#' data(goals)
#' fit0 <- gw(goals ~ offset(log(played)), data = goals)
#' summary(fit0)
#'
#' fit1 <- add1(fit0, ~ position)
#' summary(fit1)
#'
#' @export
add1.gw <- function (object, scope, test = c("none", "Chisq"), k = 2, trace = FALSE, ...) {

  safe_pchisq <- function(q, df, ...){
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }

  if (missing(scope) || is.null(scope))
    stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", scope), c("df", "AIC")))
  ans[1L, ] <- extractAIC(object, k = k, ...)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq(ns)) {
    tt <- scope[i]
    if (trace > 0) {
      cat("trying +", tt, "\n", sep = "")
      utils::flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . +", tt)), evaluate = FALSE)
    nfit <- eval(nfit, envir = env)
    ans[i + 1L, ] <- extractAIC(nfit, k = k, ...)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if (all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2L])
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- ans[, 2L] - k * ans[, 1L]
    dev <- dev[1L] - dev
    dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


#' Drop All Possible Single Terms to a GWRM Model
#'
#' Compute all the single terms in the scope argument that can be dropped from the GWRM model, fit those models and compute a table of the changes in fit.
#'
#' @param object  a fitted object of class inheriting from \code{"gw"}.
#' @param scope a formula giving the terms to be considered for dropping.
#' @param test \code{"none"}, which considers the AIC criterion, or \code{Chisq}, which is the likelihood-ratio test.
#' @param k the penalty constant in AIC / Cp.
#' @param trace	if \code{TRUE}, print out progress reports.
#' @param ...	further arguments passed to or from other methods.
#'
#' @return An object of class \code{"anova"} summarizing the differences in fit between the models.
#'
#' @importFrom stats pchisq terms drop.scope update.formula extractAIC nobs formula update as.formula
#'
#' @examples
#' data(goals)
#'
#' fit0 <- gw(goals ~ offset(log(played)), data = goals)
#' summary(fit0)
#'
#' fit1 <- step(fit0, ~ position)
#' summary(fit1)

#'
#' @export
drop1.gw<-function (object, scope, test = c("none", "Chisq"), k = 2, trace = FALSE, ...) {

  safe_pchisq <- function(q, df, ...){
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }

  tl <- attr(terms(object), "term.labels")
  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)), "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", scope), c("df", "AIC")))
  ans[1, ] <- extractAIC(object, k = k, ...)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq(ns)) {
    tt <- scope[i]
    if (trace > 0) {
      cat("trying -", tt, "\n", sep = "")
      utils::flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . -", tt)), evaluate = FALSE)
    nfit <- eval(nfit, envir = env)
    ans[i + 1, ] <- extractAIC(nfit, k = k, ...)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if (all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[1L, 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- ans[, 2L] - k * ans[, 1L]
    dev <- dev - dev[1L]
    dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

#' @export
extractAIC.gw <- function (fit, scale, k = 2, ...){
  if (fit$aic < 0 || (fit$method != "nlm" && fit$code != 0)){
    c(Inf, Inf)
  }
  else{
    n <- length(fit$residuals)
    edf <- n - fit$df.residual
    aic <- fit$aic
    c(edf, aic + (k - 2) * edf)
  }
}
