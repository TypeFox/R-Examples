# playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

## S4 experiments:
#plotCoords <- function(name, object, call, envir, ...)
#    standardGeneric("plotCoords")
#setGeneric("plotCoords", signature = c("name", "object"))
## default method:
#setMethod("plotCoords", signature(), plotCoords.default)
## dendrogram method:
#setOldClass("dendrogram")
#setMethod("plotCoords",
#          signature(name = "plot", object = "dendrogram"),
#          plotCoords.plot.dendrogram)

plotCoords <- function(name, object, call, envir, ...)
    UseMethod("plotCoords")

plotCoords.default <-
    function(name, object, call, envir, data,
             panel.args = NULL, ...)
{
    if (!is.null(panel.args)) {
        ## Lattice plot
        ## try to detect plot type, since "name" might be missing
        foo <- panel.args
        if (is.null(foo$y) && !is.null(foo$distribution)) {
            ## probably `qqmath`
            return(plotCoords.qqmath(name, object, call, envir, data,
                                     panel.args = panel.args, ...))
        }
        if (!is.null(foo$scales.3d)) {
            ## probably `cloud` (or wireframe)
            return(plotCoords.cloud(name, object, call, envir, data,
                                    panel.args = panel.args, ...))
        }
        nx <- length(foo$x)
        ny <- length(foo$y)
        if ((nx == 0) || (ny == 0))
            return(NULL)
        if (nx != ny) {
            if (nx < ny)
                foo$x <- rep(foo$x, length.out = ny)
            else foo$y <- rep(foo$y, length.out = nx)
        }
        return(foo)

    } else {
        ## non-lattice plot
        ## skip to the "plot.default" method
        plotCoords.plot.default(name, object, call, envir, data, ...)
    }
}

### BASE GRAPHICS METHODS

plotCoords.qqnorm <-
plotCoords.qqplot <-
    function(name, object, call, envir, ...)
{
    ## these return plotted coordinates in a list
    ## could also do: playDevCur()$result
    call$plot <- FALSE
    eval(call, envir)
}

plotCoords.plot <- function(name, object, call, envir, ...)
{
    ## generic plot method.
    ## dispatch on the class of 'object' (first argument to call)
    UseMethod("plotCoords.plot", object)
}

plotCoords.plot.default <- function(name, object, call, envir, data, ...)
{
    ## this is the fallback for non-lattice plots
    tmp.x <- object
    tmp.y <- eval(call[["y"]], data, envir)
    xy.coords_with_class(tmp.x, tmp.y, data = data, envir = envir)
}

plotCoords.plot.SpatialPoints <- 
plotCoords.plot.SpatialPointsDataFrame <-
    function(name, object, call, envir, data, ...)
{
    xy <- coordinates(object)
    list(x = xy[,1], y = xy[,2])
}
case.names.SpatialPoints <-
case.names.SpatialPointsDataFrame <-
    function(object, ...)
{
    case.names(as.data.frame(object), ...)
}

plotCoords.plot.dendrogram <- function(name, object, call, envir, ...)
{
    ## indices into original data for each node
    subscripts <- order.dendrogram(object)
    ## return positions in original data order
    nodePos <- order(subscripts)
    basePos <- rep(0, length(nodePos))
    horiz <- eval(call$horiz, envir)
    if (is.null(horiz)) horiz <- FALSE
    if (horiz) {
        list(x = basePos, y = nodePos)
    } else {
        list(x = nodePos, y = basePos)
    }
}
case.names.dendrogram <- function(object, ...)
{
    subscripts <- order.dendrogram(object)
    labels(object)[order(subscripts)]
}

plotCoords.plot.hclust <- function(name, object, call, envir, ...)
{
    object <- as.dendrogram(object)
    plotCoords.plot.dendrogram(name, object, call = call,
                               envir = envir, ...)
}
case.names.hclust <- function(object, ...)
    case.names(as.dendrogram(object), ...)

plotCoords.plot.mca <- function(name, object, call, envir, ...)
{
    ## only makes sense if coordinates of rows are plotted
    rows <- eval(call$rows, envir)
    if (identical(rows, FALSE)) return(NULL)
    list(x = object$rs[,1], y = object$rs[,2])
}
case.names.mca <- function(object, ...)
    case.names(object$rs)

plotCoords.plot.lm <- function(name, object, call, envir, ...)
{
    ## playwith can only handle this if it is a single plot
    ## i.e. length(which) == 1

    ## TODO
    NULL
}

plotCoords.termplot <- function(name, object, call, envir, ...)
{
    ## playwith can only handle this if it is a single plot
    ## i.e. length(terms) == 1

    ## and only makes sense if plotting residuals
    resids <- eval(call$partial.resid, envir)
    if (is.null(resids)) resids <- FALSE
    if (!isTRUE(resids)) return(NULL)

    ## TODO
    NULL
}

plotCoords.biplot <- function(name, object, call, envir, ...)
    UseMethod("plotCoords.biplot", object)

plotCoords.biplot.default <-
    function(name, object, call, envir, ...)
{
    list(x = object[,1], y = object[,2])
}

plotCoords.biplot.princomp <-
plotCoords.biplot.prcomp <- function(name, object, call, envir, ...)
{
    ## match the call to this specific method, not generic 'biplot'
    call <- match.call(stats:::biplot.prcomp, call)
    ## substitute default arguments if missing
    for (nm in c("choices", "scale", "pc.biplot")) {
        if ((nm %in% names(call)) == FALSE)
            call[nm] <- formals(stats:::biplot.prcomp)[nm]
    }
    ## extract required arguments
    x <- object
    choices <- eval(call$choices, envir)
    scale <- eval(call$scale, envir)
    pc.biplot <- eval(call$pc.biplot, envir)
    lam <- x$sdev[choices]
    ## handle different names for prcomp vs princomp:
    scores <- if (!is.null(x$scores)) x$scores else x$x
    n <- if (!is.null(x$n.obs)) x$n.obs else NROW(scores)
    lam <- lam * sqrt(n)
    if (scale != 0)
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot)
        lam <- lam/sqrt(n)
    coords <- t(t(scores[, choices])/lam)
    list(x = coords[,1], y = coords[,2])
}
case.names.prcomp <- function(object, ...)
    case.names(object$x)
case.names.princomp <- function(object, ...)
    case.names(object$scores)

### LATTICE METHODS

plotCoords.qqmath <- function(name, object, call, envir, panel.args, ...)
{
    ## based on panel.identify.qqmath
    x <- panel.args$x
    distribution <- panel.args$distribution
    groups <- panel.args$groups
    subscripts <- panel.args$subscripts
    x <- as.numeric(x)
    if (is.null(subscripts)) subscripts <- seq_along(x)
    if (!is.null(panel.args$f.value))
        return(NULL)
    distribution <-
        if (is.function(distribution)) distribution
        else if (is.character(distribution)) get(distribution)
        else eval(distribution)
    ## compute quantiles corresponding to given vector, possibly
    ## containing NA's.  The return value must correspond to the
    ## original order
    getq <- function(x)
    {
        ans <- x
        id <- !is.na(x)
        ord <- order(x[id])
        if (any(id)) ans[id] <- distribution(ppoints(sum(id)))[order(ord)]
        ans
    }
    if (is.null(groups))
    {
        return(list(x = getq(x), y = x, subscripts = subscripts))
    }
    else
    {
        allq <- rep(NA_real_, length(x))
        subg <- groups[subscripts]
        vals <- if (is.factor(groups)) levels(groups) else sort(unique(groups))
        for (i in seq_along(vals))
        {
            ok <- !is.na(subg) & (subg == vals[i])
            allq[ok] <- getq(x[ok])
        }
        return(list(x = allq, y = x, subscripts = subscripts))
    }
}

plotCoords.cloud <- function(name, object, call, envir, panel.args, ...)
{
    if (!exists("panel.identify.cloud")) return(NULL)
    idcall <- call("panel.identify.cloud", panel.args = panel.args)
    for (x in c("screen", "R.mat", "perspective", "distance", "aspect")) {
        val <- call[[x]]
        if (!is.null(val))
            idcall[x] <- list(val)
    }
    idcall$panel.3d.identify <-
        function(x, y, z, rot.mat = diag(4), distance, xlim.scaled,
                 ylim.scaled, zlim.scaled, subscripts, ...)
        {
            id <- ((x >= xlim.scaled[1]) & (x <= xlim.scaled[2]) &
                   (y >= ylim.scaled[1]) & (y <= ylim.scaled[2]) &
                   (z >= zlim.scaled[1]) & (z <= zlim.scaled[2]) &
                   !is.na(x) & !is.na(y) & !is.na(z))
            m <- ltransform3dto3d(rbind(x, y, z), rot.mat, distance)
            xpos <- m[1, ]
            ypos <- m[2, ]
            xpos[!id] <- NA
            ypos[!id] <- NA
            list(x = xpos, y = ypos, subscripts = subscripts)
        }
    eval(idcall, envir)
}

plotCoords.parallel <- function(name, object, call, envir, panel.args, ...)
{
    ## based on panel.parallel
    pcParallel <-
        function(x, y, z, subscripts, groups = NULL,
                 common.scale = FALSE,
                 lower = sapply(z, function(x) min(as.numeric(x), na.rm = TRUE)),
                 upper = sapply(z, function(x) max(as.numeric(x), na.rm = TRUE)),
                 ..., horizontal.axis = TRUE)
        {
            n.r <- ncol(z)
            n.c <- length(subscripts)
            if (is.function(lower))
                lower <- sapply(z, lower)
            if (is.function(upper))
                upper <- sapply(z, upper)
            if (common.scale) {
                lower <- min(lower)
                upper <- max(upper)
            }
            lower <- rep(lower, length = n.r)
            upper <- rep(upper, length = n.r)
            dif <- upper - lower
            if (n.r == 0) return(NULL)
            zz <- data.matrix(z[subscripts, ])
            zz <- scale(zz, center = lower, scale = dif)
            ii <- rep(factor(colnames(z), levels = colnames(z)),
                      each = n.c)
            #ii <- matrix(seq_len(n.r), byrow = TRUE,
            #             ncol = n.r, nrow = n.c)
            if (horizontal.axis) {
                list(x = zz, y = ii, subscripts = subscripts)
            } else {
                list(x = ii, y = zz, subscripts = subscripts)
            }
        }
    do.call("pcParallel", panel.args)
}

plotCoords.splom <- function(name, object, call, envir, panel.args, packet, ...)
{
    ## current viewport, restore when finished
    vp <- current.vpPath()
    on.exit({
        upViewport(0)
        if (length(vp) > 0) downViewport(vp)
    })
    upViewport(0)
    ## find the current panel and go to its viewport
    lay <- playDevCur()$tmp$currentLayout
    col <- col(lay)[lay == packet]
    row <- row(lay)[lay == packet]
    panelvp <- trellis.vpname("panel", col, row)
    downViewport(panelvp)
    ## work out plot coordinates in panelvp space
    pargs <- panel.args
    nvars <- length(pargs$z)
    ## data
    subscripts <- pargs$subscripts
    z <- data.matrix(pargs$z[subscripts, ])
    coords <- list()
    coords$x <- matrix(0, ncol = nvars * (nvars - 1),
                       nrow = length(subscripts))
    coords$y <- coords$x
    coords$subscripts <- subscripts
    #coords$x <-
    #    lapply(1:nvars, function(i)
    #       {
    #           rep(pargs$z, nvars)
    #coords$y <- rep(pargs$z, each = nvars)
    i <- 1
    for (row in 1:nvars)
        for (column in 1:nvars)
            if (row != column) {
                subpanel.name <- paste("subpanel", column, row, sep = ".")
                depth <- downViewport(subpanel.name)
                xlim <- convertX(unit(0:1, "npc"), "native", TRUE)
                ylim <- convertY(unit(0:1, "npc"), "native", TRUE)
                upViewport(depth)
                x <- z[, column]
                y <- z[, row]
                ## convert to "npc" in the subpanel
                x <- (x - xlim[1]) / (xlim[2] - xlim[1])
                y <- (y - ylim[1]) / (ylim[2] - ylim[1])
                ## convert to superpanel coordinates
                x <- (x + column - 1) + 0.5
                y <- (y + row - 1) + 0.5
                #i <- (row - 1) * nvars + column
                coords$x[,i] <- x
                coords$y[,i] <- y
                i <- i + 1
            }
    coords
}


### UTILITY

## adapted from grDevices::xy.coords
xy.coords_with_class <-
    function(x, y = NULL, recycle = TRUE, data = NULL, envir = NULL)
{
    if (is.null(y)) {
        if (is.language(x)) {
            if (inherits(x, "formula") && length(x) == 3) {
                if (!is.null(data)) {
                    y <- eval(x[[2]], data, envir)
                    x <- eval(x[[3]], data, envir)
                } else {
                    y <- eval(x[[2]], environment(x), envir)
                    x <- eval(x[[3]], environment(x), envir)
                }
            }
            else return(NULL) #stop("invalid first argument")
        }
        else if (inherits(x, "zoo")) {
            y <- if (is.matrix(x)) x[, 1] else x
            x <- stats::time(x)
        }
        else if (inherits(x, "ts")) {
            y <- if (is.matrix(x)) x[, 1] else x
            x <- stats::time(x)
        }
        else if (is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
        }
        else if (is.matrix(x) || is.data.frame(x)) {
            if (ncol(x) == 1) {
                y <- x[, 1]
                x <- seq_along(y)
            }
            else {
                y <- x[, 2]
                x <- x[, 1]
            }
        }
        else if (is.list(x)) {
            y <- x[["y"]]
            x <- x[["x"]]
        }
        else {
            y <- x
            x <- seq_along(x)
        }
    }
    if (inherits(x, "zoo")) {
        ## and y is not null
        x <- as.numeric(x)
        y <- as.numeric(y)
    }
    if (inherits(x, "POSIXt")) x <- as.POSIXct(x)
    if (length(x) != length(y)) {
        if (recycle) {
            if ((nx <- length(x)) < (ny <- length(y)))
                x <- rep(x, length.out = ny)
            else y <- rep(y, length.out = nx)
        } else stop("'x' and 'y' lengths differ")
    }
    list(x = x, y = y)
}

## MISC LABELS

case.names.zoo <-
case.names.ts <-
    function(object, ...)
{
    rep(format(stats::time(object)), length = length(object))
}

case.names.yearmon <-
case.names.yearqtr <-
case.names.Date <-
case.names.POSIXt <-
    function(object, ...)
    format(object)
