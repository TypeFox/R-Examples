#' @importFrom stats printCoefmat pnorm


summary.ah <- function(object, ...){
  options(digits =5)
    ci.lower <- object$coef - 1.96*object$se
    ci.upper <- object$coef + 1.96*object$se
    z <- object$coef/object$se
    p.value <- 2*pnorm(-abs(z))
    TAB <- cbind( coef = object$coef, se = object$se,  lower.95 = ci.lower,  upper.95= ci.upper, z = z, p.value = p.value )
  
	res <- list(call=object$call, coefficients=TAB)
	class(res) <- "summary.ah"
   
   	res
    

}

#' @export
print.summary.ah<- function(x, digits = max(getOption("digits") - 4, 4), signif.stars = getOption("show.signif.stars"),...){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
                ...)
  
}

#' @export
summary.ah.2ph <- function(object,...){
  options(digits =5)
    ci.lower <- object$coef - 1.96*object$se
    ci.upper <- object$coef + 1.96*object$se
    z <- object$coef/object$se
    p.value <- 2*pnorm(-abs(z))
    TAB <- cbind( coef = object$coef, se = object$se,  lower.95 = ci.lower,  upper.95= ci.upper, z = z, p.value = p.value)
  
	res <- list(call=object$call, coefficients=TAB)
	class(res) <- "summary.ah"
   
   	res
   

}

# print.summary.ah.2ph<- function(object, ...){
# 	cat("Call:\n")
# 	print(object$call)
# 	cat("\n")
# 	printCoefmat(object$coefficients, P.value=TRUE, has.Pvalue=TRUE, digits = 5)
# }

#' @export
print.summary.ah.2ph<- function(x, digits = max(getOption("digits") - 4, 4), signif.stars = getOption("show.signif.stars"),...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
                ...)
 
}