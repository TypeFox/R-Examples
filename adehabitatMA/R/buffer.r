"buffer" <- function(xy, x, dist)
{
    ## Verifications
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixels"))
        stop("x should inherit the class SpatialPixels")

    pfsx <- proj4string(x)
    pfsxy <- proj4string(xy)
    if (!identical(pfsx, pfsxy))
        stop("different proj4string in x and xy")

    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if (abs(gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    if (!inherits(xy, "SpatialLines")) {
        if (ncol(coordinates(xy))>2)
            stop("xy should be defined in two dimensions")
    } else {
        if (ncol(coordinates(xy)[[1]][[1]])>2)
            stop("xy should be defined in two dimensions")
    }
    meth <- "none"
    if (inherits(xy, "SpatialPoints")) {
        if (inherits(xy,"SpatialPointsDataFrame")) {
            if (ncol(xy)==1) {
                meth <- "points.id"
            } else {
                meth <- "points"
            }
        } else {
            meth <- "points"
        }
    }
    x <- as(x, "SpatialGrid")
    if (inherits(xy, "SpatialLines")) {
        meth <- "lines"
    }
    if (inherits(xy, "SpatialPolygons")) {
        meth <- "polygons"
    }


    if (meth=="none")
        stop("non convenient xy")

    if (meth == "points") {
        xy <- coordinates(xy)
        xyg <- coordinates(x)
        res <- rep(0, nrow(xyg))

        for (i in (1:nrow(xy))) {
            ii <- xy[i,]
            re <- xyg[(xyg[,1]>=(ii[1]-dist))&(xyg[,1]<=(ii[1]+dist))&(xyg[,2]>=(ii[2]-dist))&(xyg[,2]<=(ii[2]+dist)),]

            ok <- unlist(lapply(1:nrow(re), function(j) {
                k <- 0
                di <- sqrt(sum((ii-re[j,])^2))
                ok <- 0
                if (di <= dist) {
                    ok <- 1
                }
                return(ok)
            }))
            res[(xyg[,1]>=(ii[1]-dist))&(xyg[,1]<=(ii[1]+dist))&(xyg[,2]>=(ii[2]-dist))&(xyg[,2]<=(ii[2]+dist))] <- res[(xyg[,1]>=(ii[1]-dist))&(xyg[,1]<=(ii[1]+dist))&(xyg[,2]>=(ii[2]-dist))&(xyg[,2]<=(ii[2]+dist))]+ok
        }

        res <- data.frame(x=as.numeric(res>0))
        coordinates(res) <- coordinates(x)
        gridded(res) <- TRUE
        if (!is.na(pfsxy))
            proj4string(res) <- CRS(pfsxy)
        return(res)
    }
    if (meth=="points.id") {
        id <- factor(xy[[1]])
        xy <- as.data.frame(coordinates(xy))
        lixy <- split(xy, id)
        res <- do.call("data.frame", lapply(lixy, function(z) {
            bu <- buffer(SpatialPoints(z, proj4string=CRS(pfsx)), x, dist)
            return(bu[[1]])
        }))
        coordinates(res) <- coordinates(x)
        gridded(res) <- TRUE

        if (!is.na(pfsxy))
            proj4string(res) <- CRS(pfsxy)
        return(res)
    }


    if (meth=="lines") {
        gridded(x) <- TRUE
        xyg <- coordinates(x)

        li <- coordinates(xy)
        res <- apply(do.call("cbind", lapply(li, function(y) {
            apply(do.call("cbind", lapply(y, function(z) {

                gr <- gridparameters(x)

                ra <- gr[2, 2]/100
                z[,1]<-jitter(z[,1], amount=ra)
                z[,2]<-jitter(z[,2], amount=ra)

                bu <- buffer(SpatialPoints(z, proj4string=CRS(pfsx)), x, dist)


                if (nrow(gr) > 2)
                    stop("sg should be defined in two dimensions")
                if (abs(gr[1, 2] - gr[2, 2]) > get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
                    stop("the cellsize should be the same in x and y directions")
                nRow <- gr[2, 3]
                nCol <- gr[1, 3]

                carter <- matrix(0, nrow=nRow, ncol = nCol)
                xgr <- unique(xyg[,1])
                ygr <- unique(xyg[,2])

                toto <- .C("bufligr", as.double(t(z)), as.double(dist),
                           as.double(t(carter)), as.double(xgr),
                           as.double(ygr), as.integer(nCol),
                           as.integer(nRow), as.integer(nrow(z)),
                           PACKAGE="adehabitatMA")[[3]]
                toto <- c(matrix(toto, nrow = nCol, byrow=TRUE))
                toto <- toto + bu[[1]]
                ## the buffer on the line is summed to the buffer on the points
                return(toto)
            })), 1, sum)
        })), 1, sum)
        res <- data.frame(x=as.numeric(res>0))
        coordinates(res) <- coordinates(x)
        gridded(res) <- TRUE

        if (!is.na(pfsxy))
            proj4string(res) <- CRS(pfsxy)

        return(res)
    }


    if (meth=="polygons") {
        gridded(x) <- TRUE
        xyg <- coordinates(x)
        ov <- as.numeric(!is.na(over(x, xy)))
        xy <- as(xy, "SpatialLines")
        li <- coordinates(xy)
        res <- apply(do.call("cbind", lapply(li, function(y) {
            apply(do.call("cbind", lapply(y, function(z) {

                ## overlay:

                gr <- gridparameters(x)

                ra <- gr[2, 2]/100
                z[,1]<-jitter(z[,1], amount=ra)
                z[,2]<-jitter(z[,2], amount=ra)

                bu <- buffer(SpatialPoints(z, proj4string=CRS(pfsx)), x, dist)

                if (nrow(gr) > 2)
                    stop("sg should be defined in two dimensions")
                if (abs(gr[1, 2] - gr[2, 2]) > get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
                    stop("the cellsize should be the same in x and y directions")
                nRow <- gr[2, 3]
                nCol <- gr[1, 3]

                carter <- matrix(0, nrow=nRow, ncol = nCol)
                xgr <- unique(xyg[,1])
                ygr <- unique(xyg[,2])

                toto <- .C("bufligr", as.double(t(z)), as.double(dist),
                           as.double(t(carter)), as.double(xgr),
                           as.double(ygr), as.integer(nCol),
                           as.integer(nRow), as.integer(nrow(z)),
                           PACKAGE="adehabitatMA")[[3]]
                toto <- c(matrix(toto, nrow = nCol, byrow=TRUE))

                ## the buffer on the line is summed to the buffer on the points
                toto <- toto + bu[[1]]


                return(toto)
            })), 1, sum)
        })), 1, sum)
        res <- res+ov
        res <- data.frame(x=as.numeric(res>0))
        coordinates(res) <- coordinates(x)
        gridded(res) <- TRUE

        if (!is.na(pfsxy))
            proj4string(res) <- CRS(pfsxy)

        return(res)
    }
}
