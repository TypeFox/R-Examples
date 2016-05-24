plot.mgram <- function(x, pval = 0.05, xlab = "Distance", ylab = 
                "Mantel r", ...)
{
# x is the output from mgram
# pval is the p-value to be considered signficant
# ... are additional graphics parameters

        x <- x$mgram

        pval.v <- x[, 4]
        plot(x[, 1], x[, 3], type = "l", xlab = xlab, ylab = 
                ylab, ...)
        points(x[pval.v <= pval, 1], x[pval.v <= pval, 3], pch = 16)
        points(x[pval.v > pval, 1], x[pval.v > pval, 3], pch = 1)
        invisible()
}
