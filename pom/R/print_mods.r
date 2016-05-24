
print.pom <- function( x,  digits = max(3, getOption("digits") - 3), ... ){

print( x$call )
cat("Convergence Code : ", x$convergence, fill=T)
cat("\nPSI Coefficients :\n")
out <- cbind( x$psi.coefs, x$se.psi.coefs, x$psi.coefs/x$se.psi.coefs, pnorm(abs(x$psi.coefs/x$se.psi.coefs), lower.tail=F)*2 )
dimnames(out) <- list(names(x$psi.coefs), c("Estimate", "Std. Error", "Z", "Pr(>|Z|)"))
print.default(format(out, digits = digits), print.gap = 2, quote = FALSE, right = TRUE)
cat("\nP Coefficients :\n")
out <- cbind( x$p.coefs, x$se.p.coefs, x$p.coefs/x$se.p.coefs, pnorm(abs(x$p.coefs/x$se.p.coefs), lower.tail=F)*2 )
dimnames(out) <- list(names(x$p.coefs), c("Estimate", "Std. Error", "Z", "Pr(>|Z|)"))
print.default(format(out, digits = digits), print.gap = 2, quote = FALSE, right = TRUE)
cat( "\n" )
cat("AIC: ", format(x$aic, digits=digits+2), fill=TRUE)
cat("BIC: ", format(x$bic, digits=digits+2), fill=TRUE)
cat("Average PSI-hat = ", mean(x$psi.est), fill=TRUE)
cat("Average P-hat = ", mean(x$p.est), fill=TRUE)

}


print.mixed.pom <- function( x, digits = max(3, getOption("digits") - 3), ... ){

print(x$call)
cat("\n")
cat("Convergence Code : ", x$convergence, fill=T)
cat("\n")
cat("\nPSI Coefficients : \n")
out <- cbind( x$psi.coefs, x$se.psi.coefs, x$psi.coefs/x$se.psi.coefs, pnorm(abs(x$psi.coefs/x$se.psi.coefs), lower.tail=F)*2 )
dimnames(out) <- list(names(x$psi.coefs), c("Estimate", "Std. Error", "Z", "Pr(>|Z|)"))
print.default(format(out, digits = digits), print.gap = 2, quote = FALSE, right = TRUE)
cat("\n")
cat("Beta mixture parameters : \n")
print.default(format(x$p.coefs, digits = digits), print.gap = 2, quote = FALSE, right=TRUE)
cat( "\n" )
cat("AIC: ", format(x$aic, digits=digits+2), fill=TRUE)
cat("BIC: ", format(x$bic, digits=digits+2), fill=TRUE)
cat("Average PSI-hat = ", mean(x$psi.est), fill=TRUE)

}
