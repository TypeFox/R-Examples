print.sample.size.mean <-
function(x,...) {
 cat('\nsample.size.mean object: Sample size for mean estimate\n', sep="")
 if(x$call$N == Inf) cat("Without ")
 else cat("With ")
 cat("finite population correction: N=",x$call$N,", precision e=",x$call$e," and standard deviation S=",round(x$call$S,4),"\n", sep="")
 cat("\nSample size needed: ", x$n, "\n\n", sep="")
 invisible(x)
}
