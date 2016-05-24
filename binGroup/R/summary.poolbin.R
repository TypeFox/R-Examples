summary.poolbin <-
function(object, scale=object$scale, ...){
args <- list(...)
if(is.null(args$digits)) digits <- 4
else digits <- args$digits
cat("Estimation of Binomial Proportion for Pooled Data\n\n")
print(object, scale=scale, ...)
cat("\n")
cat(paste("Point estimator:",object$pt.method,"\n"))
cat(paste("CI method:",object$ci.method,"\n\n"))
cat(paste("Number of individuals:",sum(object$n * object$m),"\n"))
cat(paste("Number of pools:",sum(object$n),"\n"))
cat(paste("Number of positive pools:",sum(object$x),"\n"))
invisible(object)
}

