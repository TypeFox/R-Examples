## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer


convertToDevicePixels <-
    function(x, y)
{
    ## x and y can be unit objects or numeric
    if (!is.unit(x)) x <- unit(x, "native")
    if (!is.unit(y)) y <- unit(y, "native")
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    xy <- cbind(convertX(x, "inches", valueOnly=TRUE),
                convertY(y, "inches", valueOnly=TRUE),
                1)
    location <- xy %*% current.transform() ## inches
    dpi <- dev.size("px") / dev.size("in")
    locx <- round(dpi[1] * location[,1] * location[,3])
    locy <- round(dpi[2] * location[,2] * location[,3])
    ## convert y coordinate to have origin at top-left
    locy <- dev.size("px")[2] - locy
    list(x = locx, y = locy)
}

convertFromDevicePixels <-
    function(x.px, y.px, unitTo = "native", valueOnly = FALSE)
{
    ## convert pixels to inches
    dpi <- dev.size("px") / dev.size("in")
    x.in <- x.px / dpi[1]
    y.in <- y.px / dpi[2]
    ## convert y coordinate from origin at top-left
    y.in <- dev.size("in")[2] - y.in
    xy <- cbind(x.in, y.in, 1)
    location <- xy %*% solve(current.transform()) ## inches
    locx <- unit(location[,1] * location[,3], "inches")
    locy <- unit(location[,2] * location[,3], "inches")
    list(x = convertX(locx, unitTo, valueOnly=valueOnly),
         y = convertY(locy, unitTo, valueOnly=valueOnly))
}

inViewport <- function(x.px, y.px, viewport)
{
    ## current viewport, restore when finished
    vp <- current.vpPath()
    on.exit({
        upViewport(0)
        if (length(vp) > 0) downViewport(vp)
    })
    upViewport(0)
    if (!is.null(viewport))
        downViewport(viewport)
    ## calculate bounding box
    xy <- convertToDevicePixels(x = unit(0:1, "npc"),
                                y = unit(0:1, "npc"))
    ## viewport might be rotated, so use range of x and y
    x <- xy$x
    y <- xy$y
    ## test for point inside bounding box
    ((min(x) <= x.px) & (x.px <= max(x)) &
     (min(y) <= y.px) & (y.px <= max(y)))
}

grobBBDevicePixels <- function(grob, viewport, pad = 2)
{
    ## current viewport, restore when finished
    vp <- current.vpPath()
    on.exit({
        upViewport(0)
        if (length(vp) > 0) downViewport(vp)
    })
    upViewport(0)
    if (!is.null(viewport))
        downViewport(viewport)
    ## calculate bounding box
    if (inherits(grob, "points") ||
        inherits(grob, "lines") ||
        inherits(grob, "polyline"))
    {
        ## grobX for these refers to the convex hull,
        ## which can be bad if they are colinear
        xy <- convertToDevicePixels(x = grob$x, y = grob$y)
    } else {
        gx <- unit.c(grobX(grob, "west"),
                     grobX(grob, "east"))
        gy <- unit.c(grobY(grob, "south"),
                     grobY(grob, "north"))
        xy <- convertToDevicePixels(x = gx, y = gy)
    }
    xy$x <- range(xy$x, na.rm = TRUE) + c(-pad, pad)
    xy$y <- range(xy$y, na.rm = TRUE) + c(-pad, pad)
    xy
}

showGrobsBB <- function(...)
{
    .Deprecated("grobBoundingBoxes")
    grobBoundingBoxes(...)
}

grobBoundingBoxes <-
    function(draw = TRUE,
             gp.box = gpar(col = "yellow",
                           lwd = 5, alpha = 0.2),
             gp.text = gpar(cex = 0.75, alpha = 0.5))
{
    ## current viewport, restore when finished
    vp <- current.vpPath()
    on.exit({
        upViewport(0)
        if (length(vp) > 0) downViewport(vp)
    })
    ## get a list of all grobs in the scene
    upViewport(0)
    ## need viewports = TRUE, otherwise vpPath is empty
    objs <- as.data.frame(unclass(grid.ls(viewports = TRUE,
                                          print = FALSE)),
                          stringsAsFactors = FALSE)
    objs <- subset(objs, type == "grobListing")
    if (nrow(objs) == 0) return()
    objs$vpPath <- sub("^ROOT::", "", objs$vpPath)
    objs$vpPath[objs$vpPath == "ROOT"] <- ""
    ## calculate bounding boxes around grobs
    bblist <- list()
    length(bblist) <- nrow(objs)
    n <- 0
    for (i in seq_len(nrow(objs))) {
        vpPath <- objs$vpPath[i]
        if (vpPath == "") vpPath <- NULL
        gName <- objs$name[i]
        ## TODO: objs$gPath[i] // strict=TRUE
        grob <- grid.get(gName)
        if (inherits(grob, "nullGrob"))
            next
        ## bounding box
        bb <- grobBBDevicePixels(grob, vpPath)
        bb$name <- gName
        bb$vpPath <- vpPath
        bb$class <- class(grob)
        ## construct a display name
        displayName <- as.character(grob)
        depth <- objs$vpDepth[i] - 1
        if (depth > 0)
            displayName <- paste(paste(rep(".", depth), collapse=""),
                              displayName, sep="")
        bb$displayName <- displayName
        if (draw) {
            ## draw bounding box and grob name
            grid.rect(x=mean(bb$x), y=mean(bb$y),
                      width=diff(bb$x), height=diff(bb$y),
                      default.units="native", gp=gp.box,
                      name="TMP_BOUNDBOX")
            grid.text(label=gName, x=bb$x[1], y=bb$y[1],
                      just=c("left", "top"),
                      default.units="native", gp=gp.text,
                      name="TMP_BOUNDBOX")
        }
        n <- n + 1
        bblist[[n]] <- bb
    }
    length(bblist) <- n
    ## remove annotations from display list
    ## but do not redraw it yet (so still visible)
    if (draw)
        grid.remove("TMP_BOUNDBOX", global=TRUE, redraw=FALSE)
    invisible(bblist)
}

identifyGrob <- function(xy.pixels = grid.locator(), classes = NULL)
{
    ## bounding boxes for all grobs in the scene
    bblist <- grobBoundingBoxes(draw = FALSE)
    if (length(bblist) == 0) stop("No grobs found.")
    force(xy.pixels)
    if (is.null(xy.pixels)) return(NULL)
    xy.pixels <- lapply(xy.pixels, as.numeric)
    x.px <- xy.pixels$x
    y.px <- xy.pixels$y
    grobNames <- NULL
    for (i in length(bblist):1) {
        obj <- bblist[[i]]
        x <- obj$x
        y <- obj$y
        if (!is.null(classes)) {
            if (!any(classes %in% obj$class))
                next
        }
        if ((min(x) <= x.px) && (x.px <= max(x)) &&
            (min(y) <= y.px) && (y.px <= max(y))) {
            grobNames <- c(grobNames, obj$name)
        }
    }
    grobNames
}
