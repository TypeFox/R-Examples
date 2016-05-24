print.logistftest <-
function(x, ...){
# x ... object of class logistftest
 print(x$call)
 cat("Model fitted by", x$method, "\n\nFactors fixed as follows:\n")
 print(x$testcov)
 LL <- 2 * diff(x$loglik)
 out <- c(x$loglik[1], x$loglik[2], LL/2)
 names(out) <- c("Restricted model", "Full model", "difference")
 cat("\nLikelihoods:\n")
 print(out)
 cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", x$prob, "\n", sep = "")
 invisible(x)
}

