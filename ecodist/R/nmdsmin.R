nmds.min <- function(x, dims=0)

{
# returns the minimum-stress configuration from nmds output.
# (Results from nmds)
# if dims==0, returns the overall lowest-stress configuration 
# Otherwise, returns the lowest-stress configuration of dimensionality dims

    if(dims == 0) {
        x.min <- x$conf[x$stress == min(x$stress)]
    }
    else {
        x.dims <- sapply(x$conf, ncol)
        x$conf <- x$conf[x.dims == dims]
        x$stress <- x$stress[x.dims == dims]
        x$r2 <- x$r2[x.dims == dims]
        x.min <- x$conf[x$stress == min(x$stress)]
    }
    cat("Minimum stress for given dimensionality: ", x$stress[which.min(x$stress)], "\n")
    cat("r^2 for minimum stress configuration: ", x$r2[which.min(x$stress)], "\n")
    x.min <- x.min[[1]]
    x.min <- data.frame(x.min)
    attr(x.min, "stress") <- x$stress[which.min(x$stress)]
    attr(x.min, "r2") <- x$r2[which.min(x$stress)]
    x.min
}
