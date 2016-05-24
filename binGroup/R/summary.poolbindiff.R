summary.poolbindiff <-
function(object, scale=object$scale, ...){
args <- list(...)
if(is.null(args$digits)) digits <- 4
else digits <- args$digits
cat("Estimation of Difference of Binomial Proportions for Pooled Data\n\n")
print(object, scale=scale, ...)
cat("\n")
cat(paste("Point estimator:",object$pt.method,"\n"))
cat(paste("CI method:",object$ci.method,"\n\n"))
p1 <- pooledBin(object$x1,object$m1,object$n1,scale=object$scale)
p2 <- pooledBin(object$x2,object$m2,object$n2,scale=object$scale)
summat <- matrix(c(scale*p1$p, scale*p1$lcl, scale*p1$ucl, scale, sum(object$n1 * object$m1), sum(object$n1), sum(object$x1),
             scale*p2$p, scale*p2$lcl, scale*p2$ucl, scale, sum(object$n2 * object$m2), sum(object$n2), sum(object$x2)), nrow=2,ncol=7, byrow=TRUE)
summat <- round(summat, digits=digits)          
dimnames(summat) <- list(c("Population 1","Population 2"), c("PointEst","Lower","Upper","Scale","Individuals","Pools","Positive Pools"))         
if(scale == 1) summat <- summat[,-4]
print(summat)
invisible(object)
}

