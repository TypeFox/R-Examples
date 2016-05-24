#' @title Model comparison test for fitted cdfqr models
#' @description Likelihood Ratio Tests for fitted cdfqr Objects.
#' @aliases anova.cdfqr
#' @param object The fitted cdfqr model.
#' @param ... One or more cdfqr model objects for model comparison.
#' @param test The model comparision test, currently only 'LRT' is implemented.
#' @examples
#' data(cdfqrExampleData)
#' fit_null <- cdfquantreg(crc99 ~ 1 | 1, 't2','t2', data = JurorData)
#' fit_mod1 <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' anova(fit_null, fit_mod1)
#' 
#' @method anova cdfqr
#' @export

anova.cdfqr <- function(object,..., test = "LRT") {
  models <- list(...)
  if (length(models)) 
    object <- c(list(object), models)
  else warning("Sorry, the current test is for model comparison only. Please supply more than one models for comparison.")
  
  nmodels <- length(object)
  object <- object[order(-as.numeric(lapply(object, function(x) x$df.residual)))]
  
  if (nmodels == 1) 
    warning("Sorry, the single model anova test is not available yet. Please run null model, and compare the null model with the present model")
  
  resdf <- as.numeric(lapply(object, function(x) x$df.residual))
  resdev <- as.numeric(lapply(object, function(x) -2*x$logLik))
  
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(formula(x)),  collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "-2Loglik", 
                                        "Df", "Deviance"))
  
  table <- stats::stat.anova(table = table, test = "LRT", scale = 1, 
                      df.scale = min(resdf), n = max(resdf))
  
  names(table) <-  c("Resid. Df", "-2Loglik",  "Df", "LR stat","Pr(>Chi)")
  class(table) <- c("anova", "data.frame")
  attr(table, "heading") <- c("Likelihood ratio tests \n")
  
  table
}