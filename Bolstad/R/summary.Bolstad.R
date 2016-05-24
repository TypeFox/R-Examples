#' Summarizing Bayesian Multiple Linear Regression
#' 
#' \code{summary} method for output of \code{\link{bayes.lm}}.
#' 
#' 
#' @param object an object of "\code{Bolstad}" that is the result of a call to \code{\link{bayes.lm}}
#' @param \dots any further arguments to be passed to \code{print}
#' 
#' 
#' @seealso The function to fit the model \code{\link{bayes.lm}}
#' 
#' The function \code{\link{coef}} to extract the matrix of posterior means along with standard errors and t-statistics.
#' 
#' @export

summary.Bolstad = function(object, ...) {
  if(length(class(object)) == 2 && all(grepl("Bolstad|lm", class(object)))){
    
    getTermLabels = function(x){ attr(x$terms, "term.labels") }
    
    z = object
    
    ans = list(rank = z$rank,
               call = z$call,
               terms = c("(Intercept)", getTermLabels(z)),
               coef = z$coefficients,
               std.err = diag(z$post.var),
               prior = z$prior,
               residuals = as.vector(z$residuals),
               res.df = z$df.residual, ...)
    class(ans) = "summary.Bolstad"
    ans
  }else{
    summary.default(object, ...)
  } 
}

