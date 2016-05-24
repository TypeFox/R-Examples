## Print function for sscopu objects
print.sscopu <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## dimension
    cat("Dimemsion: ",dim(x$basis)[2],".",sep="")
    cat("\n\n")
    order <- x$order
    if (is.null(order)) order <- 2
    cat("Maximum order of interaction: ",order,".",sep="")
    cat("\n\n")
    if (x$symmetry) {
        cat("The fit is symmetric, invariant to variable permutation.")
        cat("\n\n")
    }
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for sshzd2d objects
print.sshzd2d <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    if (x$symmetry) {
        cat("The fit is symmetric, with a common marginal hazard and a symmetric copula.")
        cat("\n\n")
    }
    ## terms
    cat("Terms in hzd1:\n")
    print.default(c(x$hzd1$terms$labels,x$hzd1$lab.p))
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n")
    print.default(x$hzd1$desc)
    cat("\n")
    cat("Terms in hzd2:\n")
    print.default(c(x$hzd2$terms$labels,x$hzd2$lab.p))
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n")
    print.default(x$hzd2$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}
