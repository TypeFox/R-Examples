findmaxasc <- function(asc)
{
    if (!inherits(asc,"asc"))
        stop("asc should be of class \"asc\"")
    xyc <- getXYcoords(asc)
    asc2 <- asc
    asc[is.na(asc)] <- -9999

    toto <- .C("findmaxgrid",as.double(t(asc)), as.integer(nrow(asc)),
               as.integer(ncol(asc)), PACKAGE="adehabitat")[[1]]

    toto <- c(matrix(toto, ncol=ncol(asc), byrow=TRUE))
    xy <- expand.grid(xyc$x,xyc$y)
    xyb <- xy[which(toto>0.5),]
    xyb <- xyb[!is.na(join.asc(xyb,asc2)),]
    return(xyb)
}

