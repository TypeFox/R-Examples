mimage <- function(x,  var = names(slot(x,"data")),
                   col = gray((240:1)/256), mfrow = NULL)

{
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("x should be of class SpatialPixelsDataFrame")

    slot(x, "data") <- slot(x, "data")[,var]
    if (is.null(mfrow))
        mfrow <- n2mfrow(length(var))
    na <- names(slot(x,"data"))
    opar <- par(mfrow=mfrow, mar=c(0,0,2,0))
    on.exit(par(opar))
    lapply(na, function(i) {
        image(x, i, col=col)
        title(main=i)
        box()
    })
    invisible(NULL)
}
