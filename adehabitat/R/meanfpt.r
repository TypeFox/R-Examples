"meanfpt" <- function(f, graph=TRUE)
{
    ## Verifications
    if (!inherits(f, "fipati"))
        stop("x should be of class 'fipati'")

    ## Graphical setting
    if (graph)
        opar <- par(mfrow=n2mfrow(length(f)))

    ## gets the radii
    s <- attr(f, "radii")

    ## Computes and plots the mean FPT for each radius and for each burst
    soso <- lapply(f, function(y) {
        so <- apply(y,2,function(z) mean(z, na.rm=TRUE))
        if (graph)
            plot(s, so, ty="l", xlab="scale", ylab="Mean of FPT",
                 main=attr(y,"burst"))
        return(so)
    })

    ## output
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

