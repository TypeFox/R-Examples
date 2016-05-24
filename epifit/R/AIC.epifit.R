##' Function for calculating Akaike's \sQuote{An Information Criterion} (AIC) from epifit object.
##'
##' Function called from generic function AIC in \pkg{stats} when the argument is epifit object.
##' @param object a fitted epifit object.
##' @param ... not used in this version, only for compatibility purpose with generic function \code{AIC} currently.
##' @param k numeric, the \emph{penalty} per parameter to be used; the default \code{k = 2} is the classical AIC.
##' @return a numeric AIC value.
##' @seealso \code{\link[stats]{AIC}}.
##' @examples
##' library(survival)
##' # The simplest test data set from coxph function
##' test1 <- list(time=c(4,3,1,1,2,2,3),
##'               status=c(1,1,1,0,1,1,0),
##'               x=c(0,2,1,1,1,0,0),
##'               sex=c(0,0,0,0,1,1,1))
##' AIC(coxph(Surv(time,status)~x + strata(sex), data=test1))
##' modelexpr <- "cox(time,status)/strata(sex)~exp(beta*x)"
##' AIC(epifit(modelexpr=modelexpr, data=test1))
##' @export
AIC.epifit <- function(object, ..., k = 2){
  object$loglik[2]*(-2) + length(object$coefficients)*k
}
