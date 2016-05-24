##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

tileplot <-
    function(x, data = NULL, aspect = "iso",
             prepanel = "prepanel.default.xyplot",
             panel = "panel.voronoi", ...)
{
    foo <- levelplot(x, data = data, aspect = aspect,
                     panel = panel, prepanel = prepanel, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

## panel function to draw Voronoi mosaic
panel.voronoi <-
    function(x, y, z, subscripts = TRUE, at = pretty(z),
             points = TRUE, border = "transparent",
             na.rm = FALSE, win.expand = 0.07, use.tripack = FALSE,
             ...,
             col.regions = regions$col, alpha.regions = regions$alpha)
{
    ## We need either tripack (better? but weird license) or
    ## deldir. Go with deldir unless explicitly requested.
    if (use.tripack) {
        if (!requireNamespace("tripack", quietly = TRUE))
            stop("The 'use.tripack=TRUE' option requires the 'tripack' package to be installed.")
    } else {
        if (!requireNamespace("deldir", quietly = TRUE))
            stop("This function requires the 'deldir' package to be installed.")
    }
    ## find subset of points to use
    x0 <- x[subscripts]
    y0 <- y[subscripts]
    z0 <- z[subscripts]
    ## throw away NAs, but keep originals for panel.xyplot()
    ok <- complete.cases(x0, y0)
    if (na.rm) ok <- ok & !is.na(z0)
    x <- x0[ok]
    y <- y0[ok]
    z <- z0[ok]
    if (!any(is.finite(z))) return()
    ## strip duplicated locations, with warning
    dup <- duplicated(cbind(x, y))
    if (any(dup)) {
        warning(paste("Ignoring", sum(dup),
                      "cases of duplicated locations"))
        x <- x[!dup]
        y <- y[!dup]
        z <- z[!dup]
    }
    ## compute bounds
    data.rg <- list(x = extendrange(x, f = win.expand),
                    y = extendrange(y, f = win.expand))
    bounds <- c(data.rg$x, data.rg$y)
    #panel.rg <- lapply(current.panel.limits(), sort)
    #bounds <- c(max(panel.rg$x[1], data.rg$x[1]),
    #            min(panel.rg$x[2], data.rg$x[2]),
    #            max(panel.rg$y[1], data.rg$y[1]),
    #            min(panel.rg$y[2], data.rg$y[2]))
    ## check if any points in visible plot region
    #if (is.unsorted(bounds[1:2]))
    #    bounds[1:2] <- panel.rg$x
    #if (is.unsorted(bounds[3:4]))
    #    bounds[3:4] <- panel.rg$y
    if (use.tripack) {
        xy <- data.frame(x = x, y = y)
        ## add dummy points to ensure that voronoi polygons are finite
        dummies <- data.frame(x = c(-1,-1,1,1), y = c(-1,1,-1,1)) * 10 * max(abs(xy))
        xy <- rbind(xy, dummies)
        tiles <- tripack::voronoi.polygons(tripack::voronoi.mosaic(xy, duplicate = "error"))
    } else {
        ## NB: the 'rw' argument as subset of data is bad because
        ## need to take corresponding subset of z !
        ## (but not easy to work out what that is)

        #set <- ((bounds[1] < x) & (x < bounds[2]) &
        #        (bounds[3] < y) & (y < bounds[4]))
        #x <- x[set]
        #y <- y[set]
        #z <- z[set]
        tiles <- deldir::tile.list(deldir::deldir(x, y, rw = bounds))
        tiles <- lapply(tiles, function(p) as.data.frame(p[c("x", "y")]))
    }
    ## draw it as one composite polygon
    polydata <- do.call("rbind", tiles)
    regions <- trellis.par.get("regions")
    zcol <- level.colors(z, at, col.regions, colors = TRUE)
    grid.polygon(polydata[,1], polydata[,2],
                 id.lengths = sapply(tiles, nrow),
                 default.units = "native",
                 gp = gpar(fill = zcol, col = border,
                 alpha = alpha.regions))
    if (points) {
        panel.xyplot(x0, y0, ...)
    }
}

panel.levelplot.points <-
    function(x, y, z, subscripts = TRUE, at = pretty(z),
             shrink, labels, label.style, contour, region, ## (all ignored)
             pch = 21, col.symbol = "#00000044",
             ...,
             col.regions = regions$col,
             fill = NULL) ## (ignored)
{
    regions <- trellis.par.get("regions")
    zcol <- level.colors(z, at, col.regions, colors = TRUE)
    x <- x[subscripts]
    y <- y[subscripts]
    zcol <- zcol[subscripts]
    ## panel.xyplot does the work (can handle 'type' argument, etc)
    panel.xyplot(x, y, fill = zcol, pch = pch,
                 col.symbol = col.symbol, ...)
}
