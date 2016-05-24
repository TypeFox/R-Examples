.charhull <- function(xy)
{
    xy <- as.data.frame(xy)
    names(xy) <- c("x","y")
    tri <- deldir(xy[,1], xy[,2])
    tri <- triang.list(tri)
    tri <- lapply(1:length(tri), function(i) {
        rbind(tri[[i]], tri[[i]][1,])
    })
    area <- sapply(tri, function(x) {
        .arcpspdf(SpatialPolygons(list(Polygons(list(Polygon(x[,2:3])), 1))))
    })
    tri <- tri[order(area)]
    area <- area[order(area)]
    df <- data.frame(area=area, percent=100*c(c(1:length(area))/length(area)))
    sp <- SpatialPolygons(lapply(1:length(tri), function(i) {
        Polygons(list(Polygon(as.matrix(tri[[i]][,2:3]))), i)
    }))
    res <- SpatialPolygonsDataFrame(sp, df)
    return(res)
}

CharHull <-  function(xy, unin = c("m", "km"), unout = c("ha", "m2", "km2"),
              duplicates = c("random", "remove"), amount = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class SpatialPoints")
    pfs <- proj4string(xy)
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy) != 1) {
            warning("xy contains more than one column. id ignored")
            m <- 1
            id <- rep("a", nrow(coordinates(xy)))
        }
        else {
            id <- xy[[1]]
            m <- 2
        }
    }  else {
        m <- 1
        id <- rep("a", nrow(coordinates(xy)))
    }
    id <- factor(id)
    duplicates <- match.arg(duplicates)
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    xy <- as.data.frame(coordinates(xy))
    lixy <- split(xy, id)
    res <- list()

    for (i in 1:length(lixy)) {
        x <- as.matrix(lixy[[i]])
        levv <- factor(apply(x, 1, paste, collapse = " "))
        if (duplicates == "remove") {
            x <- as.data.frame(x)
            x <- as.matrix(do.call("rbind", lapply(split(x, levv),
                                                   function(z) z[1, ])))
        }
        else {
            x <- as.data.frame(x)
            sam1 <- x[, 1] - jitter(x[, 1], 1, amount)
            sam2 <- x[, 1] - jitter(x[, 1], 1, amount)
            lixyt <- split(x, levv)
            lisam1 <- split(sam1, levv)
            lisam2 <- split(sam2, levv)
            x <- as.matrix(do.call("rbind", lapply(1:length(lixyt),
                function(a) {
                  z <- lixyt[[a]]
                  if (nrow(z) > 1) {
                    z[, 1] <- z[, 1] + lisam1[[a]]
                    z[, 2] <- z[, 2] + lisam2[[a]]
                  }
                  return(z)
                })))
        }


        resa <- .charhull(x)

        pol <- as(resa, "SpatialPolygons")
        lip <- list(pol[1])
        for (j in 2:nrow(resa)) {
            poo <- rbind(pol[j], lip[[j-1]])
            pls <- slot(poo, "polygons")
            pls1 <- lapply(pls, maptools::checkPolygonsHoles)
            slot(poo, "polygons") <- pls1
            lip[[j]] <- rgeos::gUnionCascaded(poo, id = rep(j,
                                                   length(row.names(poo))))
        }
        are <- .arcpspdf(lip[[1]])
        for (j in 2:length(lip)) {
            are[j] <- .arcpspdf(lip[[j]])
        }
        spP <- do.call("rbind", lip)
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
        df <- data.frame(area = are, percent = resa[[2]])
        resa <- SpatialPolygonsDataFrame(spP, df)
        if (!is.na(pfs))
            proj4string(resa) <- CRS(pfs)
        res[[i]] <- resa
    }


    if (m == 1) {
        return(res[[1]])
    }
    else {
        names(res) <- names(lixy)
        class(res) <- "MCHu"
        return(res)
    }
}
