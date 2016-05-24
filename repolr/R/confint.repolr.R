confint.repolr <-
function (object, parm, level = 0.95, robust.var = TRUE, ...) {
 cf <- object$coefficients
 pnames <- names(cf)
 if (missing(parm) == TRUE){ 
   parm <- pnames
 } else if (is.numeric(parm) == TRUE){ 
   parm <- pnames[parm]
 }
 a <- (1 - level)/2
 a <- c(a, 1 - a)
 fac <- qnorm(a)
 pct <- paste(format(100 * a, trim = TRUE, 
          scientific = FALSE, digits = 3),"%")
 ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
 if(robust.var == TRUE){
  ses <- sqrt(diag(object$robust.var))[parm]
 } else {
  ses <- sqrt(diag(object$naive.var))[parm]
 }
 ci[] <- cf[parm] + ses %o% fac
 ci
}
