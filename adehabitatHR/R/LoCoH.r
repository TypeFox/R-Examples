
## k-Locoh

LoCoH.k <- function(xy, k=5, unin = c("m", "km"),
                    unout = c("ha", "m2", "km2"),
                    duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    pfs <- proj4string(xy)
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")

    duplicates <- match.arg(duplicates)
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    m <- 1
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        xy <- coordinates(xy)
        if (k > nrow(xy))
            stop("too large number of neighbours")

        ## Management of duplicates
        levv <- factor(apply(xy,1,paste, collapse=" "))
        if (duplicates=="remove") {
            xy <- as.data.frame(xy)
            xy <- as.matrix(do.call("rbind", lapply(split(xy, levv), function(x) x[1,])))
        } else {
            xy <- as.data.frame(xy)
            sam1 <- xy[,1] - jitter(xy[,1], 1, amount)
            sam2 <- xy[,1] - jitter(xy[,1], 1, amount)
            lixy <- split(xy, levv)
            lisam1 <- split(sam1, levv)
            lisam2 <- split(sam2, levv)
            xy <- as.matrix(do.call("rbind",
                                    lapply(1:length(lixy),
                                           function(i){
                                               x <- lixy[[i]]
                                               if (nrow(x)>1) {
                                                   x[,1] <- x[,1]+lisam1[[i]]
                                                   x[,2] <- x[,2]+lisam2[[i]]
                                               }
                                               return(x)
                                           })))
        }


        ind <- 1:nrow(xy)

        ## identification of the clusters
        oo <- do.call("rbind",lapply(1:nrow(xy), function(i) {
            di <- unlist(lapply(1:nrow(xy), function(j) {
                sqrt(sum((xy[j,]-xy[i,])^2))
            }))
            ind2 <- ind[order(di)]
            return(ind2[1:k])
        }))

        ## identification of the coordinates of the MCP
        pol <- lapply(1:nrow(xy), function(i) {
            ff <- do.call("rbind", lapply(1:k, function(j) xy[oo[i,j],]))
            pol <- ff[chull(ff[,1], ff[,2]),]
            return(pol)
        })

        ## computes area of the polygons
        ar <- unlist(lapply(1:nrow(xy), function(i) {
            .arcpspdf(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]], pol[[i]][1,]))), i))))
        }))

        ## sort everything according to the area:
        xy <- xy[order(ar),]
        pol <- pol[order(ar)]
        oo <- oo[order(ar),]
        ar <- ar[order(ar)]

        ## then, "incremental" union:
        lip <- list(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[1]], pol[[1]][1,]))), 1))))
        n <- k
        dej <- unlist(oo[1,])

        for (i in 2:nrow(xy)) {
            poo <- rbind(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]],pol[[i]][1,]))), i))), lip[[i-1]])
            pls <- slot(poo, "polygons")
            pls1 <- lapply(pls, maptools::checkPolygonsHoles)
            slot(poo, "polygons") <- pls1
            lip[[i]] <- rgeos::gUnionCascaded(poo, id=rep(i, length(row.names(poo))))
            dej <- c(dej, unlist(oo[i,!(unlist(oo[i,])%in%dej)]))
            n[i] <- length(dej)
        }

        ## Compute the area
        are <- .arcpspdf(lip[[1]])
        for (i in 2:nrow(xy)) {
            are[i] <- .arcpspdf(lip[[i]])
        }

        ## And the results, as a SpatialPolygonDataFrame object
        spP <- do.call("rbind", lip)

        ## The data frame:
        if (unin == "m") {
            if (unout == "ha")
                are <- are/10000
            if (unout == "km2")
                are <- are/1e+06
        }
        if (unin == "km") {
            if (unout == "ha")
                are <- are * 100
            if (unout == "m2")
                are <- are * 1e+06
        }

        df <- data.frame(area=are, percent=100*n/nrow(xy))
        res <- SpatialPolygonsDataFrame(spP, df)
        if (!is.na(pfs))
            proj4string(res) <- CRS(pfs)
        return(res)
    } else {
        id <- as.data.frame(xy)[[1]]
        xy <- as.data.frame(coordinates(xy))
        lixy <- split(xy, id)
        res <- lapply(lixy, function(x) {
            LoCoH.k(SpatialPoints(x, proj4string=CRS(as.character(pfs))),
                    k, unin, unout)
        })
        class(res) <- "MCHu"
        return(res)
    }
}


LoCoH.k.area <- function(xy, krange, percent=100, unin = c("m", "km"),
                         unout = c("ha", "m2", "km2"),
                         duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    if (any(percent>100)) {
	warning("The MCP is estimated using all relocations (percent>100)")
	percent<-100
    }

    duplicates <- match.arg(duplicates)
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    m <- 1
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        re <- unlist(lapply(krange, function(k) {
            re <- LoCoH.k(xy, k, unin, unout,
                          duplicates, amount)
            df <- as.data.frame(re)
            ind <- 1:nrow(df)
            i <- min(ind[df$percent>=percent])
            return(df$area[i])
        }))
        res <- data.frame(k=krange,area=re)
        plot(res, ty="b", xlab="k", ylab="area")
        invisible(re)
    } else {
        re <- as.data.frame(do.call("rbind", lapply(krange, function(k) {
            re <- LoCoH.k(xy, k, unin,
                          unout, duplicates, amount)
            unlist(lapply(re, function(kk) {
                df <- as.data.frame(kk)
                ind <- 1:nrow(df)
                i <- min(ind[df$percent>=percent])
                return(df$area[i])
            }))
        })))
        opar <- par(mfrow=n2mfrow(ncol(re)))
        lapply(1:ncol(re), function(i) {
            plot(krange, re[,i], ty="b", xlab="k", ylab="area", main=names(re)[i])

        })
        par(opar)
        invisible(re)
    }
}



## Rasterization

MCHu.rast <- function(x, spdf, percent=100)
{
    if ((!inherits(x, "SpatialPolygonsDataFrame"))&(!inherits(x, "MCHu")))
        stop("x should have been generated by LoCoH.* or clusthr")
    if (!inherits(spdf, "SpatialPixels"))
        stop("spdf should inherit the class SpatialPixels")
    pfs2 <- proj4string(spdf)

    gridded(spdf) <- TRUE
    gr <- gridparameters(spdf)
    if (nrow(gr) > 2)
        stop("spdf should be defined in two dimensions")


    if (inherits(x, "SpatialPolygonsDataFrame")) {
        pfs1 <- proj4string(x)
        if (!identical(pfs1,pfs2))
            stop("x and spdf do not have the same proj4string")
        tmp <- rep(0, nrow(as.data.frame(spdf)))
        spdf <- as(spdf,"SpatialPixels")
        df <- as.data.frame(x)
        if (max(df$percent)<(max(percent)))
            stop("The number of isolated points was too large to allow the computation of the specified isopleth\nTry to decrease percent")
        ind <- 1:nrow(df)
        lb <- unlist(lapply(percent, function(p) {
            min(ind[df$percent>=p])
        }))
        x <- x[lb,]

        for (i in nrow(as.data.frame(x)):1) {
            uu <- !is.na(over(spdf,geometry(x[i,])))
            tmp[uu] <- df$percent[i]
        }
        tmp[tmp==0] <- NA
        tmp2 <- data.frame(tmp)
        coordinates(tmp2) <- coordinates(spdf)
        gridded(tmp2) <- TRUE
        if (!is.na(pfs1))
            proj4string(tmp2) <- CRS(pfs1)
        return(tmp2)
    } else {
        re <- lapply(x, function(y) {
            as.data.frame(MCHu.rast(y, spdf, percent))[,1]
        })
        re <- do.call("data.frame",re)
        coordinates(re) <- coordinates(spdf)
        gridded(re) <- TRUE
        if (!is.na(pfs2))
            proj4string(re) <- CRS(pfs2)
        return(re)
    }
}





### R LoCoH

LoCoH.r <- function(xy, r, unin = c("m", "km"),
                    unout = c("ha", "m2", "km2"),
                    duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    pfs <- proj4string(xy)
    m <- 1
    duplicates <- match.arg(duplicates)
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        xy <- coordinates(xy)
        ind <- 1:nrow(xy)

        ## Management of duplicates
        levv <- factor(apply(xy,1,paste, collapse=" "))
        if (duplicates=="remove") {
            xy <- as.data.frame(xy)
            xy <- as.matrix(do.call("rbind", lapply(split(xy, levv), function(x) x[1,])))
        } else {
            xy <- as.data.frame(xy)
            sam1 <- xy[,1] - jitter(xy[,1], 1, amount)
            sam2 <- xy[,1] - jitter(xy[,1], 1, amount)
            lixy <- split(xy, levv)
            lisam1 <- split(sam1, levv)
            lisam2 <- split(sam2, levv)
            xy <- as.matrix(do.call("rbind",
                                    lapply(1:length(lixy),
                                           function(i){
                                               x <- lixy[[i]]
                                               if (nrow(x)>1) {
                                                   x[,1] <- x[,1]+lisam1[[i]]
                                                   x[,2] <- x[,2]+lisam2[[i]]
                                               }
                                               return(x)
                                           })))
        }


        ## identification of the clusters
        oo <- lapply(1:nrow(xy), function(i) {
            di <- unlist(lapply(1:nrow(xy), function(j) {
                sqrt(sum((xy[j,]-xy[i,])^2))
            }))
            ind2 <- ind[order(di)]
            di <- di[order(di)]
            return(ind2[di <= r])
        })

        ## identification of the coordinates of the MCP
        pol <- lapply(1:nrow(xy), function(i) {
            ff <- do.call("rbind",
                          lapply(1:length(oo[[i]]), function(j) xy[oo[[i]][j],]))
            pol <- ff[chull(ff[,1], ff[,2]),]
            if (is.null(nrow(pol)))
                pol <- matrix(pol, nrow=1)
            return(pol)
        })

        ## computes area of the polygons
        ar <- unlist(lapply(1:nrow(xy), function(i) {
            if (nrow(pol[[i]])>2) {
                return(.arcpspdf(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]], pol[[i]][1,]))), i)))))
            } else {
                return(0)
            }
        }))

        ## Computes the number of relocations in each polygon
        num <- unlist(lapply(oo, length))

        ## sort everything according to the number of relocations,
        ## and then area
        ind <- order(-num, ar)
        xy <- xy[ind,]
        pol <- pol[ind]
        ar <- ar[ind]
        oo <- oo[ind]
        num <- num[ind]
        lone <- length(ar[ar<.Machine$double.eps])
        max <- nrow(xy)-lone
        if (max==0)
            stop("the distance is too small: there were only isolated points\nPlease increase the value of r")

        ## then, "incremental" union:
        lip <- list(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[1]], pol[[1]][1,]))), 1))))
        dej <- oo[[1]]
        n <- num[1]
        for (i in 2:nrow(xy)) {
            if (nrow(pol[[i]])>2) {
                poo <- rbind(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]],pol[[i]][1,]))), i))), lip[[i-1]])
            } else {
                poo <- lip[[i-1]]
            }
            pls <- slot(poo, "polygons")
            pls1 <- lapply(pls, maptools::checkPolygonsHoles)
            slot(poo, "polygons") <- pls1

            lip[[i]] <- rgeos::gUnionCascaded(poo, id=rep(i, length(row.names(poo))))
            dej <- c(dej, oo[[i]][!(oo[[i]]%in%dej)])
            n[i] <- length(dej)
        }

        ## Compute the area
        are <- .arcpspdf(lip[[1]])
        for (i in 2:nrow(xy)) {
            are[i] <- .arcpspdf(lip[[i]])
        }

        ## And the results, as a SpatialPolygonDataFrame object
        spP <- do.call("rbind", lip)

        ## The data frame:
        n[n>max] <- max
        if (unin == "m") {
            if (unout == "ha")
                are <- are/10000
            if (unout == "km2")
                are <- are/1e+06
        }
        if (unin == "km") {
            if (unout == "ha")
                are <- are * 100
            if (unout == "m2")
                are <- are * 1e+06
        }


        df <- data.frame(area=are, percent=100*n/nrow(xy))
        res <- SpatialPolygonsDataFrame(spP, df)
        if (!is.na(pfs))
            proj4string(res) <- CRS(pfs)
        return(res)

    } else {
        id <- as.data.frame(xy)[,1]
        xy <- as.data.frame(coordinates(xy))
        lixy <- split(xy, id)
        res <- lapply(lixy, function(x) {
            LoCoH.r(SpatialPoints(x, proj4string=CRS(as.character(pfs))),
                    r, unin, unout)
        })
        class(res) <- "MCHu"
        return(res)
    }
}




LoCoH.r.area <- function(xy, rrange, percent=100, unin = c("m", "km"),
                         unout = c("ha", "m2", "km2"),
                         duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    duplicates <- match.arg(duplicates)
    m <- 1
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        re <- unlist(lapply(rrange, function(k) {
            re <- LoCoH.r(xy, k, unin,
                          unout,
                          duplicates, amount)
            df <- as.data.frame(re)
            if (max(df$percent)< percent)
                stop("The number of isolated points was too large to allow the computation of the specified isopleth\nTry to decrease percent")
            ind <- 1:nrow(df)
            i <- min(ind[df$percent>=percent])
            return(df$area[i])
        }))
        res <- data.frame(k=rrange,area=re)
        plot(res, ty="b", xlab="r", ylab="area")
        invisible(re)
    } else {
        re <- as.data.frame(do.call("rbind", lapply(rrange, function(k) {
            re <- LoCoH.r(xy, k, unin,
                          unout,
                          duplicates, amount)
            return(unlist(lapply(re, function(kk) {
                df <- as.data.frame(kk)
                if (max(df$percent)< percent)
                    stop("The number of isolated points was too large to allow the computation of the specified isopleth\nTry to decrease percent or to increase the min of ranger")
                ind <- 1:nrow(df)
                i <- min(ind[df$percent>=percent])
                return(df$area[i])
            })))
        })))
        opar <- par(mfrow=n2mfrow(ncol(re)))
        lapply(1:ncol(re), function(i) {
            plot(rrange, re[,i], ty="b", xlab="r", ylab="area", main=names(re)[i])

        })
        par(opar)
        invisible(re)
    }
}








### A LoCoH

LoCoH.a <- function(xy, a, unin = c("m", "km"),
                    unout = c("ha", "m2", "km2"),
                    duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    pfs <- proj4string(xy)
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    duplicates <- match.arg(duplicates)
    m <- 1
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        xy <- coordinates(xy)
        ind <- 1:nrow(xy)


        ## Management of duplicates
        levv <- factor(apply(xy,1,paste, collapse=" "))
        if (duplicates=="remove") {
            xy <- as.data.frame(xy)
            xy <- as.matrix(do.call("rbind",
                                    lapply(split(xy, levv), function(x) x[1,])))
        } else {
            xy <- as.data.frame(xy)
            sam1 <- xy[,1] - jitter(xy[,1], 1, amount)
            sam2 <- xy[,1] - jitter(xy[,1], 1, amount)
            lixy <- split(xy, levv)
            lisam1 <- split(sam1, levv)
            lisam2 <- split(sam2, levv)
            xy <- as.matrix(do.call("rbind",
                                    lapply(1:length(lixy),
                                           function(i){
                                               x <- lixy[[i]]
                                               if (nrow(x)>1) {
                                                   x[,1] <- x[,1]+lisam1[[i]]
                                                   x[,2] <- x[,2]+lisam2[[i]]
                                               }
                                               return(x)
                                           })))
        }


        ## identification of the clusters
        oo <- lapply(1:nrow(xy), function(i) {
            di <- unlist(lapply(1:nrow(xy), function(j) {
                sqrt(sum((xy[j,]-xy[i,])^2))
            }))
            ind2 <- ind[order(di)]
            di <- cumsum(di[order(di)])
            return(ind2[di <= a])
        })

        ## identification of the coordinates of the MCP
        pol <- lapply(1:nrow(xy), function(i) {
            ff <- do.call("rbind",
                          lapply(1:length(oo[[i]]), function(j) xy[oo[[i]][j],]))
            pol <- ff[chull(ff[,1], ff[,2]),]
            if (is.null(nrow(pol)))
                pol <- matrix(pol, nrow=1)
            return(pol)
        })

        ## computes area of the polygons
        ar <- unlist(lapply(1:nrow(xy), function(i) {
            if (nrow(pol[[i]])>2) {
                return(.arcpspdf(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]], pol[[i]][1,]))), i)))))
            } else  {
                return(0)
            }
        }))

        ## Computes the number of relocations in each polygon
        num <- unlist(lapply(oo, length))

        ## sort everything according to the number of relocations,
        ## and then area
        ind <- order(-num, ar)
        xy <- xy[ind,]
        pol <- pol[ind]
        ar <- ar[ind]
        oo <- oo[ind]
        num <- num[ind]

        ## then, "incremental" union:
        lip <- list(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[1]], pol[[1]][1,]))), 1))))
        dej <- oo[[1]]
        n <- num[1]
        for (i in 2:nrow(xy)) {
            if (nrow(pol[[i]])>3) {
                poo <- rbind(SpatialPolygons(list(Polygons(list(Polygon(rbind(pol[[i]],pol[[i]][1,]))), i))), lip[[i-1]])
            } else {
                poo <- lip[[i-1]]
            }
            pls <- slot(poo, "polygons")
            pls1 <- lapply(pls, maptools::checkPolygonsHoles)
            slot(poo, "polygons") <- pls1

            lip[[i]] <- rgeos::gUnionCascaded(poo, id=rep(i, length(row.names(poo))))
            dej <- c(dej, oo[[i]][!(oo[[i]]%in%dej)])
            n[i] <- length(dej)
        }

        ## Compute the area
        are <- .arcpspdf(lip[[1]])
        for (i in 2:nrow(xy)) {
            are[i] <- .arcpspdf(lip[[i]])
        }

        ## And the results, as a SpatialPolygonDataFrame object
        spP <- do.call("rbind", lip)

        ## The data frame:
        df <- data.frame(area=are, percent=100*n/nrow(xy))
        if (unin == "m") {
            if (unout == "ha")
                are <- are/10000
            if (unout == "km2")
                are <- are/1e+06
        }
        if (unin == "km") {
            if (unout == "ha")
                are <- are * 100
            if (unout == "m2")
                are <- are * 1e+06
        }


        res <- SpatialPolygonsDataFrame(spP, df)
        if (!is.na(pfs))
            proj4string(res) <- CRS(pfs)
        return(res)

    } else {
        id <- as.data.frame(xy)[,1]
        xy <- as.data.frame(coordinates(xy))
        lixy <- split(xy, id)
        res <- lapply(lixy, function(x) {
            LoCoH.a(SpatialPoints(x, proj4string=CRS(as.character(pfs))),
                    a, unin, unout)
        })
        class(res) <- "MCHu"
        return(res)
    }
}


LoCoH.a.area <- function(xy, arange, percent=100, unin = c("m", "km"),
                         unout = c("ha", "m2", "km2"),
                         duplicates=c("random","remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class \"SpatialPoints\"")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    duplicates <- match.arg(duplicates)
    m <- 1
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)==1) {
            m <- 2
        } else {
            warning("xy contains more than one column; no id considered")
        }
    }

    if (m==1) {
        re <- unlist(lapply(arange, function(k) {
            re <- LoCoH.a(xy, k, unin,
                          unout,
                          duplicates, amount)
            df <- as.data.frame(re)
            if (max(df$percent)< percent)
                stop("The number of isolated points was too large to allow the computation of the specified isopleth\nTry to decrease percent")
            ind <- 1:nrow(df)
            i <- min(ind[df$percent>=percent])
            return(df$area[i])
        }))
        res <- data.frame(k=arange,area=re)
        plot(res, ty="b", xlab="a", ylab="area")
        invisible(re)
    } else {
        re <- as.data.frame(do.call("rbind", lapply(arange, function(k) {
            re <- LoCoH.a(xy, k, unin,
                          unout,
                          duplicates, amount)
            return(unlist(lapply(re, function(kk) {
                df <- as.data.frame(kk)
                if (max(df$percent)< percent)
                    stop("The number of isolated points was too large to allow the computation of the specified isopleth\nTry to decrease percent or to increase the min of ranger")
                ind <- 1:nrow(df)
                i <- min(ind[df$percent>=percent])
                return(df$area[i])
            })))
        })))
        opar <- par(mfrow=n2mfrow(ncol(re)))
        lapply(1:ncol(re), function(i) {
            plot(arange, re[,i], ty="b", xlab="a", ylab="area", main=names(re)[i])

        })
        par(opar)
        invisible(re)
    }
}









print.MCHu <- function(x, ...)
{
    if (!inherits(x,"MCHu"))
        stop("x should have been generated by the function LoCoH.* or clusthr")
    cat("********** Multiple convex hull Home range of several Animals ************\n\n")
    cat("This object is a list with one component per animal.\n")
    cat("Each component is an object of class SpatialPolygonsDataFrame\n")
    cat("The home range has been estimated for the following animals:\n")
    print(names(x))
}


plot.MCHu <- function(x, percent="all", points=NULL, ...)
{
    if ((!inherits(x, "MCHu"))&(!inherits(x, "SpatialPolygonsDataFrame")))
        stop("x should have been generated by LoCoH.* or clusthr")

    if (inherits(x, "SpatialPolygonsDataFrame")) {

        if (!is.null(points)) {
            if (!inherits(points, "SpatialPoints"))
                stop("points should inherit the class SpatialPoints")
        }

        df <- as.data.frame(x)
        if (percent[1]!="all") {
            percent <- sort(percent)
            ind <- 1:nrow(df)
            lb <- unlist(lapply(percent, function(p) {
                min(ind[df$percent>=p])
            }))
            x <- x[lb,]
            df <- df[lb,]
        }
        co <- grey(c(1:nrow(df))/nrow(df))
        plot(x, col=co, ...)
        if (!is.null(points))
            plot(points, add=TRUE)

    } else {

        if (!is.null(points)) {
            if (!inherits(points, "SpatialPointsDataFrame")) {
                stop("points should be of class SpatialPointsDataFrame")
            }
            po <- factor(as.data.frame(points)[,1])
            if (!all(names(x)%in%levels(po)))
                stop("Not all home ranges have an associated id in points")
        }

        opar <- par(mfrow=n2mfrow(length(x)))
        on.exit(par(opar))
        tmp <- lapply(1:length(x), function(u) {



            opar <- par(mar=c(0.1,0.1,2,0.1))
            plot.MCHu(x[[u]], percent, ...)
            title(names(x)[u])
            box()
            if (!is.null(points))
                plot(points[po==names(x)[u],], add=TRUE)
            par(opar)
        })
    }
}
