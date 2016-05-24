#' The confint function adapted for vcovHC
#' 
#' The confint.lm uses the t-distribution as the default
#' confidence interval estimator. When there is reason to believe that 
#' the normal distribution is violated an alternative approach using
#' the \code{\link[sandwich]{vcovHC}} may be more suitable.
#'   
#' @param object The regression model object, either an ols or lm object 
#' @param HC_type See options for \code{\link[sandwich]{vcovHC}}
#' @param t_distribution A boolean for if the t-distribution should be used
#'  or not. Defaults to FALSE. According to Cribari-Nieto and Lima's study from
#'  2009 this should not be the case.
#' @param ... Additional parameters that are passed on to 
#'  \code{\link[sandwich]{vcovHC}}
#' @return \code{matrix} A matrix (or vector) with columns giving lower and 
#'  upper confidence limits for each parameter. These will be labelled as 
#'  (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' 
#' @example inst/examples/confint_robust_example.R
#' 
#' @inheritParams stats::confint
#' @references \href{http://www.tandfonline.com/doi/abs/10.1080/00949650801935327}{
#'  F. Cribari-Neto and M. da G. A. Lima, 
#'  "Heteroskedasticity-consistent interval estimators", 
#'  Journal of Statistical Computation and Simulation, 
#'  vol. 79, no. 6, pp. 787-803, 2009.}
#' @export
confint_robust <- function(object, parm, level = 0.95, 
    HC_type="HC3", t_distribution = FALSE,...){
  cf <- coef(object); pnames <- names(cf)
  if(missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  
  a <- (1-level)/2; a <- c(a, 1-a)
  pct <- paste(format(100 * a, 
                      trim = TRUE, 
                      scientific = FALSE, 
                      digits = 3), 
               "%")
  if (t_distribution)
    fac <- qt(a, object$df.residual)
  else
    fac <- qnorm(a)
  ci <- array(NA, 
      dim = c(length(parm), 2L), 
      dimnames = list(parm, pct))
  ses <- sqrt(diag(vcovHC(object, type=HC_type, ...)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}
