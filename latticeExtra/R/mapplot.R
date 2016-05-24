
prepanel.mapplot <- function(x, y, map, ...) {
    list(xlim = range(map$x, finite = TRUE),
         ylim = range(map$y, finite = TRUE))
}

panel.mapplot <-
    function(x, y, map, breaks, colramp,
             exact = FALSE, lwd = 0.5, ...)
{
    names(x) <- tolower(as.character(y))
    mapnames <- tolower(map$names)
    mapval <- x[mapnames]
    xmatched <- names(x) %in% mapnames
    if (any(!xmatched) && !exact) {
        ## lump sub-regions together (strip name after ':')
        mapnames <- gsub(":.*$", "", mapnames)
        ## only replace values which did not match exactly
        mapval <- ifelse(is.na(mapval), x[mapnames], mapval)
        xmatched <- xmatched | (names(x) %in% mapnames)
    }
    if (any(!xmatched))
        warning(sum(!xmatched), " unmatched regions: ",
                toString(y[!xmatched], width = 60))
    interval <-
        cut(mapval, breaks = breaks,
            labels = FALSE, include.lowest = TRUE)
    col.regions <- colramp(length(breaks) - 1)
    col <- col.regions[interval]
    panel.polygon(map, col = col, lwd = lwd, ...)
}


mapplot <- function(x, data, ...) UseMethod("mapplot")


mapplot.formula <-
    function(x, data, map, outer = TRUE,
             prepanel = prepanel.mapplot,
             panel = panel.mapplot,
             aspect = "iso",
             legend = NULL,
             breaks, cuts = 30,
             colramp = colorRampPalette(brewer.pal(n = 11, name = "Spectral")),
             colorkey = TRUE,
             ## col.regions,
             ## alpha.regions,
             ...)
{
    colrampNULL <- is.null(colramp)
    if (is.null(colramp))
        colramp <- function(n)
            colorRampPalette(trellis.par.get("regions")$col)(n)
    ccall <- match.call()
    ccall$data <- data
    ccall$map <- map
    ccall$outer <- outer
    ccall$prepanel <- prepanel
    ccall$panel <- panel
    ccall$aspect <- aspect
    ccall$legend <- legend
    ccall$colramp <- colramp
    ccall$default.scales <- list(x = list(tck = 1), y = list(tck = 1))
    ccall[[1]] <- quote(lattice::dotplot)
    ans <- eval(ccall, parent.frame())
    if (missing(breaks))
    {
        x <- unlist(lapply(ans$panel.args, "[[", "x"))
        breaks <-
            if (is.factor(x)) seq_len(1 + nlevels(x)) - 0.5
            else do.breaks(range(x, finite = TRUE), cuts)
    }
##     regions <- trellis.par.get("col.regions")
##     if (missing(col.regions)) col.regions <- regions$col
##     if (missing(alpha.regions)) alpha.regions <- regions$alpha
    if (colorkey) {
        keydef <- list(at = breaks, draw = FALSE)
        if (!colrampNULL) keydef$col <- colramp(length(breaks))
        ans <-
            update(ans,
                   breaks = breaks,
                   legend = lattice:::updateList(ans$legend,
                                                 list(right = 
                                                      list(fun = draw.colorkey,
                                                           args = list(key = keydef)))))
    }
    ans$call <- sys.call(sys.parent())
    ans$call[[1]] <- quote(mapplot)
    ans
}

