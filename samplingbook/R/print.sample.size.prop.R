print.sample.size.prop <-
function(x,...) {
 cat('\nsample.size.prop object: Sample size for proportion estimate\n', sep="")
 if(x$call$N == Inf) cat("Without ")
 else cat("With ")
 cat("finite population correction: N=",x$call$N,", precision e=",x$call$e," and expected proportion P=",round(x$call$P,4),"\n", sep="")
 cat("\nSample size needed: ", x$n, "\n\n", sep="")
 invisible(x)
}
