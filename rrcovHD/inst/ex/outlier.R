## from mvoutlier
.pbb <- function(map = "bss.background", add.plot = FALSE, ...)
{
    all = get(map)
    xrange = c(min(all[, 1], na.rm = TRUE), max(all[, 1], na.rm = TRUE))
    yrange = c(min(all[, 2], na.rm = TRUE), max(all[, 2], na.rm = TRUE))
    if (!add.plot) {
        plot(1, 1, xlim = xrange, ylim = yrange, xlab = "", ylab = "",  ...)
    }
    lines(all, ...)
}


    bsstop <- NULL
##  VT::24.11.2013 - not possible to read files from the data directory - they are compressed
##  - use data() instead
##
##    load(system.file("data","bsstop.rda", package="mvoutlier"))
##    load(system.file("data","bss.background.rda", package="mvoutlier"))

    library(mvoutlier)
    data(bsstop)
    data(bss.background)

    x=bsstop[,5:14]

    ## visualize multivariate outliers in the map
    op <- par(mfrow=c(2,2))

    ## identify multivariate outliers
    x.out=OutlierPCOut(x)
    cat("\nPCOUT: ", length(which(!getFlag(x.out))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.out) + 2)
    title("PCOUT")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))


    ## compare with outlier detection based on MCD:
    x.mcd=OutlierMahdist(x, control="mcd")
    cat("\nMCD: ", length(which(!getFlag(x.mcd))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO,bsstop$YCOO,pch=16,col=getFlag(x.mcd) + 2)
    title("MCD")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    ## compare with outlier detection based on PCDIST:
    x.pcdist <- OutlierPCDist(x, explvar=0.99)
    cat("\nPCDIST: ", length(which(!getFlag(x.pcdist))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.pcdist)+2)
    title("PCDIST")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    ## compare with outlier detection based on SIGN1:
    x.sign1 <- OutlierSign1(x)
    cat("SIGN1: ", length(which(!getFlag(x.sign1))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.sign1)+2)
    title("SIGN1")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    par(op)
