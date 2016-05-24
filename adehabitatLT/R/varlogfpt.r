"varlogfpt" <- function(f, graph=TRUE)
{
    ## Verifications
    if (!inherits(f, "fipati"))
        stop("f should be of class 'fipati'")

    ## graphical settings
    if (graph)
        opar <- par(mfrow=n2mfrow(length(f)))

    ## The radii
    s <- attr(f, "radii")

    ## Computation oof the variance for each radius and each burst
    soso <- lapply(f, function(y) {
        so <- apply(y,2,function(z) var(log(z), na.rm=TRUE))
        if (graph)
            plot(s, so, ty="l", xlab="scale", ylab="Variance of log(FPT)",
                 main=attr(y,"burst"))
        return(so)
    })

    ## Output
    soso <- as.data.frame(do.call("rbind",soso))
    row.names(soso) <- unlist(lapply(f, function(z) attr(z, "burst")))
    names(soso) <- paste("r",1:ncol(soso), sep="")
    attr(soso, "radii") <- attr(f,"radii")
    if (graph)
        par(opar)
    if (graph) {
        invisible(soso)
    } else {
        return(soso)
    }
}

