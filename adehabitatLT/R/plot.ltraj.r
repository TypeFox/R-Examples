plot.ltraj <- function (x, id = unique(unlist(lapply(x, attr, which="id"))),
                        burst = unlist(lapply(x, attr, which="burst")),
                        spixdf = NULL, spoldf = NULL, xlim = NULL,
                        ylim = NULL, colspixdf = gray((240:1)/256),
                        colspoldf = "green", addpoints = TRUE,
                        addlines = TRUE, perani = TRUE, final = TRUE, ...)
{
    if (!is.null(spoldf)) {
        if (!inherits(spoldf, "SpatialPolygons"))
            stop("spoldf should inherit the class SpatialPolygons")
    }
    if (!is.null(spixdf)) {
        if (!inherits(spixdf, "SpatialPixelsDataFrame"))
            stop("spixdf should inherit the class SpatialPixelsDataFrame")
    }

    if (!inherits(x, "ltraj"))
        stop("x should be an object of class ltraj")

    ## Select id and burst
    x <- x[id=id]
    x <- x[burst=burst]
    typeII <- attr(x, "typeII")

    ## Remove NA
    x <- lapply(x, function(i) {
        jj <- i[!is.na(i$x),]
        attr(jj, "id") <- attr(i,"id")
        attr(jj, "burst") <- attr(i,"burst")
        return(jj)
    })
    class(x) <- c("ltraj","list")
    attr(x, "typeII") <- typeII
    attr(x,"regular") <- is.regular(x)

    ## keeps only the coordinates
    uu <- lapply(x, function(i) {
        i[,c("x","y")]
    })

    ## Which kind of plot?
    if (!perani)
        idc <- "burst"
    else idc <- "id"

    id <- unique(unlist(lapply(x, function(i) {
        attr(i, idc)
    })))

    if (length(id)>1)
        opar <- par(mar = c(0.1, 0.1, 2, 0.1), mfrow = n2mfrow(length(id)))

    ## xlim is the limit of the graph
    if (is.null(xlim)) {
        if (perani) {
            idtt <- unique(id(x))
            oo <- lapply(idtt,
                         function(i) unlist(lapply(x[id=i], function(j) j$x)))
            maxxl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
            xlim <- lapply(oo, function(ki) c(min(ki), min(ki)+maxxl))
        } else {
            maxxl <- max(unlist(lapply(x,
                                       function(ki) diff(range(ki$x)))))
            xlim <- lapply(x, function(ki) c(min(ki$x), min(ki$x)+maxxl))
        }
    }  else {
        xlim <- split(rep(xlim, length(id)), gl(length(id),2))
    }

    if (is.null(ylim)) {
        if (perani) {
            idtt <- unique(id(x))
            oo <- lapply(idtt,
                         function(i) unlist(lapply(x[id=i], function(j) j$y)))
            maxyl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
            ylim <- lapply(oo, function(ki) c(min(ki), min(ki)+maxyl))
        } else {
            maxyl <- max(unlist(lapply(x, function(ki) diff(range(ki$y)))))
            ylim <- lapply(x, function(ki) c(min(ki$y), min(ki$y)+maxyl))
        }
    } else {
        ylim <- split(rep(ylim, length(id)), gl(length(id),2))
    }
    names(xlim) <- id
    names(ylim) <- id

    for (i in id) {
        if (!is.null(spixdf)) {
            if (length(id)==1) {
                image(spixdf, col = colspixdf, xlim = xlim[i][[1]],
                      ylim = ylim[i][[1]])
            } else {
                image(spixdf, col = colspixdf, xlim = xlim[i][[1]],
                      ylim = ylim[i][[1]])
                title(main = i)
            }
        } else {
            if (length(id)==1) {
                plot(1, 1, type = "n", asp = 1,
                     xlim = xlim[i][[1]],
                     ylim = ylim[i][[1]], ...)
            } else {
                plot(1, 1, type = "n", asp = 1,
                     xlim = xlim[i][[1]],
                     ylim = ylim[i][[1]], axes = (length(id)==1),
                     main = i, ...)
            }
        }
        if (length(id)>1)
            box()
        if (!is.null(spoldf)) {
            plot(spoldf, add=TRUE, col = colspoldf)
        }
        if (addlines) {
            if (idc=="burst") {
                lines(x[burst=i][[1]]$x, x[burst=i][[1]]$y)
            } else {
                xtmp <- x[id=i]
                for (j in 1:length(xtmp)) {
                    lines(xtmp[[j]]$x, xtmp[[j]]$y)
                }
            }
        }
        if (addpoints) {
            if (idc=="burst") {
                points(x[burst=i][[1]]$x, x[burst=i][[1]]$y,
                       pch = 21, col = "black", bg = "white")
            } else {
                xtmp <- x[id=i]
                for (j in 1:length(xtmp)) {
                    points(xtmp[[j]]$x, xtmp[[j]]$y,
                           pch = 21, col = "black", bg = "white")
                }
            }
        }
        if (final) {
            if (idc=="burst") {
                points(x[burst=i][[1]]$x[c(1, length(x[burst=i][[1]]$x))],
                       x[burst=i][[1]]$y[c(1, length(x[burst=i][[1]]$y))],
                       pch = c(2,14), col = c("blue","red"), cex=2, lwd=2)
            } else {
                xtmp <- x[id=i]
                for (j in 1:length(xtmp)) {
                    points(xtmp[[j]]$x[c(1, length(xtmp[[j]]$x))],
                           xtmp[[j]]$y[c(1, length(xtmp[[j]]$x))],
                           pch = c(2,14), col = c("blue","red"), cex=2, lwd=2)
                }
            }

        }
    }
    if (length(id)>1)
        par(opar)
}




