"mcp.area" <- function(xy, percent = seq(20,100, by=5),
                       unin=c("m", "km"),
                       unout=c("ha", "km2", "m2"), plotit = TRUE)
{
    ## Verifications
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    id <- factor(as.data.frame(xy)[,1])

    if (any(percent>100)) {
	stop("Incorrect value for percent")
    }


    tmp <- lapply(percent, function(x) {
        mcp(xy, percent=x, unin, unout)
    })

    ar <- as.data.frame(do.call("rbind", lapply(tmp, function(x) {
        as.data.frame(x)$area
    })))
    names(ar) <- as.data.frame(tmp[[1]])$id


    ## output
    row.names(ar)<-percent
    class(ar)<-c("hrsize", "data.frame")
    if (plotit)
        plot(ar)
    return(ar)
}
