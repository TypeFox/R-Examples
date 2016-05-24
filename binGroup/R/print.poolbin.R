print.poolbin <-
function(x, scale=x$scale, ...){
args <- list(...)
if(is.null(args$digits)) digits <- 4
else digits <- args$digits
p <- round(scale*x$p,digits)
lcl <- round(scale*x$lcl, digits)
ucl <- round(scale*x$ucl, digits)
mat <- matrix(c(p,lcl,ucl,scale),nrow=1) # really to match Hmisc's binconf()
dimnames(mat) <- list(c(""),c("PointEst","Lower","Upper","Scale"))
if(scale == 1) mat <- mat[,-4]
print(mat)
invisible(x)
}

