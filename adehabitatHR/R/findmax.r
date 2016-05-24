findmax <- function(x)
{
    if (!inherits(x,"SpatialPixelsDataFrame"))
        stop("x should be of class \"SpatialPixelsDataFrame\"")
    gridded(x) <- TRUE
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")

    x2 <- x
    x <- as.image.SpatialGridDataFrame(x[,1])
    z <- x$z
    z[is.na(z)] <- -9999

    toto <- .C("findmaxgrid",as.double(t(z)), as.integer(nrow(z)),
               as.integer(ncol(z)), PACKAGE="adehabitatHR")[[1]]

    toto <- c(matrix(toto, ncol=ncol(z), byrow=TRUE))
    xy <- expand.grid(x$x,x$y)
    xyb <- xy[which(toto>0.5),]
    ov <- z[over(SpatialPoints(xyb, proj4string=CRS(proj4string(x2))), geometry(x2))]
    xyb <- SpatialPoints(xyb[!is.na(ov),])
    return(xyb)
}

