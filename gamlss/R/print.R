##---------------------------------------------------------------------------------------
##                    print.gamlss 
##---------------------------------------------------------------------------------------
print.gamlss <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{  
   # digits <-  max(3, getOption("digits") - 3) 
    cat("\nFamily: ", deparse(x$family), "\n") 
    cat("Fitting method:", deparse(x$method), "\n") 
   # cat("\nCall: ", deparse(x$call), "\n\n")
    cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
    cat("Mu Coefficients")
    if (is.character(co <- x$contrasts)) 
        cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
            paste, collapse = "="), "]")
    cat(":\n")
    if ("mu"%in%x$parameters) 
       {  
       print.default(format(coef(x,"mu"), digits = digits), print.gap = 2, quote = FALSE)
       }
    if ("sigma"%in%x$parameters) 
       {
       cat("Sigma Coefficients:\n")
       print.default(format(coef(x,"sigma"), digits = digits),print.gap = 2,quote = FALSE)
       }
    if ("nu"%in%x$parameters) 
       {  
       cat("Nu Coefficients:\n")
       print.default(format(coef(x,"nu"), digits = digits),print.gap = 2, quote = FALSE)
       }
    if ("tau"%in%x$parameters) 
       { 
       cat("Tau Coefficients:\n")
       print.default(format(coef(x,"tau"), digits = digits), print.gap = 2, quote = FALSE)
       }
    cat("\n Degrees of Freedom for the fit:", x$df.fit, "Residual Deg. of Freedom  ", 
        x$df.residual, "\n")
    cat("Global Deviance:    ", format(signif(x$G.deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), 
        "\n            SBC:    ", format(signif(x$sbc)), "\n")
    invisible(x)
}
#---------------------------------------------------------------------------------------
