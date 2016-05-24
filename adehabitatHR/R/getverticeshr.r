getverticeshr.estUD <- function(x, percent=95, ida=NULL,
                                unin = c("m", "km"),
                                unout = c("ha", "km2", "m2"),
                                standardize=FALSE, ...)
{
    ## Verifications
    if (!inherits(x,"estUD"))
        stop("x should be of class \"estUD\"")
    if (inherits(x,"estUD")&(is.null(ida)))
        ida <- "homerange"
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    pfs <- proj4string(x)
    if (!slot(x, "vol"))
        x<-getvolumeUD(x, standardize=standardize)

    ## Check that all the contour are within the study area limits
    tmp <- x[[1]]
    gp <- gridparameters(x)[,3]
    tmpm <- matrix(tmp, ncol=gp[2], nrow=gp[1], byrow = TRUE)
    ma <- min(c(tmpm[c(1:nrow(tmpm)), c(1,ncol(tmpm))],
                tmpm[c(1,nrow(tmpm)), c(1:ncol(tmpm))]))
    if (any(percent>=ma))
        stop(paste("The grid is too small to allow the estimation of home-range.\nYou should rerun kernelUD with a larger extent parameter", sep=""))


    if (length(percent)>1)
        stop("percent should be of length 1")

    xyma <- coordinates(x)
    xyl <- list(x=unique(xyma[,1]), y=unique(xyma[,2]))
    ud <- as.image.SpatialGridDataFrame(x[,1])$z

    re <- contourLines(x = xyl$x,
                       y = xyl$y,
                       ud, nlevels = 1,
                       levels = percent)

    ## identify whether there is a hole in the list:
    areaa <- unlist(lapply(re, function(y) {
        ttmp <- cbind(y$x,y$y)
        ttmp <- rbind(ttmp, ttmp[1,])
        .arcpspdf(SpatialPolygons(list(Polygons(list(Polygon(ttmp)), 1))))
    }))

    spatpol <- do.call("cbind",lapply(1:length(re), function(i) {
        y <- re[[i]]
        zz <- cbind(y$x,y$y)
        zz <- rbind(zz,zz[1,])
        tmp <- SpatialPolygons(list(Polygons(list(Polygon(zz)),
                                             as.character(i))), proj4string=CRS(pfs))
        return(!is.na(over(x, tmp)))
    }))
    spatpol <- as.data.frame(spatpol)
    hol <- unlist(lapply(1:ncol(spatpol), function(i) {
        all(apply(data.frame(spatpol[spatpol[,i],-i]), 1, any))
    }))
    areaa <- sum(areaa*sign(as.numeric(!hol)-0.5))

    ii <- SpatialPolygons(list(Polygons(lapply(1:length(re), function(i) {
        y <- re[[i]]
        zz <- cbind(y$x,y$y)
        zz <- rbind(zz,zz[1,])
        return((Polygon(zz, hole=hol[i])))
    }), ida)))

    if (unin == "m") {
        if (unout == "ha")
            areaa <- areaa/10000
        if (unout == "km2")
            areaa <- areaa/1e+06
    }
    if (unin == "km") {
        if (unout == "ha")
            areaa <- areaa * 100
        if (unout == "m2")
            areaa <- areaa * 1e+06
    }

    dff <- data.frame(id=ida, area=areaa)
    row.names(dff) <- ida
    ii <- SpatialPolygonsDataFrame(ii, dff)
    if (!is.na(pfs))
        proj4string(ii) <- CRS(pfs)
    return(ii)
}



getverticeshr.estUDm <- function(x, percent=95, whi=names(x),
                                 unin = c("m", "km"),
                                 unout = c("ha", "km2", "m2"),
                                 standardize=FALSE, ...)
{
    if (!inherits(x,"estUDm"))
        stop("x should be of class \"estUDm\"")
    x <- x[whi]

    unin <- match.arg(unin)
    unout <- match.arg(unout)

    res <- do.call("rbind", lapply(1:length(x), function(i) {
        getverticeshr(x[[i]], percent, ida=names(x)[i], unin,
                      unout, standardize, ...)
    }))
    return(res)
}




getverticeshr.MCHu <- function(x, percent=95, whi=names(x), ...)
{
    ## Verifications
    if (!inherits(x, "MCHu"))
        stop("x should be of class \"MCHu\"")
    x <- x[whi]

    ## gets the home range for a given level (percent)
    res2 <- lapply(1:length(x), function(r) {
        y <- x[[r]]
        y <- spChFIDs(y, paste(names(x)[r], 1:nrow(y), sep="."))
        per <- slot(y,"data")$percent
        dif <- percent-per
        cons <- c(1:length(per))[dif>0]
        if (length(cons)!=0) {
            whi <- max(cons)
        } else {
            stop(paste(percent,"% contour could not be created.\n More data points are probably needed."))
        }
        return(y[whi,])
    })

    res2 <- do.call("rbind", res2)
    res2 <- spChFIDs(res2, names(x))
    if (!is.na(proj4string(x[[1]])))
        proj4string(res2) <- CRS(proj4string(x[[1]]))
    df <- slot(res2, "data")
    slot(res2, "data") <- data.frame(id = row.names(df), area=df$area)
    return(res2)
}


getverticeshr <- function(x, percent=95, ...)
{
    UseMethod("getverticeshr")
}

getverticeshr.default <- function(x, percent=95, ...)
{
    stop("no method provided for this class")
}


