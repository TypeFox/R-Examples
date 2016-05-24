##
## Copyright (c) 2010 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

horizonplot <- function(x, data, ...)
    UseMethod("horizonplot")

horizonplot.default <-
    function(x, data = NULL, ...,
             nbands = 3L,
             horizonscale = NA,
             origin = function(y) na.omit(y)[1],
             colorkey = FALSE, legend = NULL,
             panel = panel.horizonplot,
             prepanel = prepanel.horizonplot,
             col.regions = brewer.pal(n = 2 * nbands, name = "RdYlBu"),
             strip = FALSE, strip.left = TRUE,
             par.strip.text = list(cex = 0.6),
             colorkey.digits = 3,
             layout = c(1, NA),
             groups = NULL,
             default.scales =
               list(y = list(relation = "free", axs = "i",
                             draw = FALSE, tick.number = 2)))
{
    if (!is.null(groups))
        stop("'groups' does not work in this plot")
    ans <- xyplot(x, data = data, ...,
                  origin = origin, horizonscale = horizonscale,
                  panel = panel, prepanel = prepanel,
                  col.regions = col.regions,
                  strip = strip, strip.left = strip.left,
                  par.strip.text = par.strip.text,
                  layout = layout,
                  default.scales = default.scales,
                  nbands = nbands)
    ans$call <- match.call()
    ## add colorkey
    if (isTRUE(colorkey)) colorkey <- list()
    if (is.list(colorkey))
    {
        bands.at <- seq(-nbands, nbands)
        if (ans$y.scales$relation == "same") {
            origin <- ans$y.limits[1]
            horizonscale <- diff(ans$y.limits)
        }
        if (is.na(horizonscale)) {
            ## labels <- expression(
            ##    - 3 * Delta[i], - 2 * Delta[i], - 1 * Delta[i], 0,
            ##    + 1 * Delta[i], + 2 * Delta[i], + 3 * Delta[i], 0)
            labels <- parse(text = sprintf("%+d * Delta[i]", bands.at))
            labels[nbands + 1] <- if (is.numeric(origin)) origin else "origin"
        }
        else {
            if (is.numeric(origin)) {
                labels <- round(origin + bands.at * horizonscale, colorkey.digits)
            } else {
                labels <- sprintf("%+g", round(bands.at * horizonscale, colorkey.digits))
                labels[nbands + 1] <- "origin"
            }
        }
        ii <- round(seq(1, length(col.regions), length.out = 2 * nbands))
        colorkey <-
            modifyList(list(col = col.regions[ii], at = bands.at,
                            labels = list(labels = labels, at = bands.at)),
                       colorkey)
        space <- colorkey$space
        if (is.null(space)) space <- "right"
        if (is.null(legend)) legend <- list()
        legend[[space]] <- list(fun = "draw.colorkey",
                                args = list(colorkey)) 
        ans <- update(ans, legend = legend)
    }
    ans
}


panel.horizonplot <-
    function(x, y, ..., border = NA,
             nbands = 3L,
             col.regions = brewer.pal(n = 2 * nbands, name = "RdYlBu"),
             origin) ## catch origin, don't pass to panel.xyarea!
{
    origin <- current.panel.limits()$ylim[1]
    scale <- diff(current.panel.limits()$ylim)
    ## ordered for drawing, from least extreme to most extreme
    #sections <- c(0, -1, 1, -2, 2, -3) ## these are the lower bounds
    sections <- as.vector(rbind(seq_len(nbands)-1, -seq_len(nbands)))
    #ii <- round(((sections + 3) / 5) * (length(col.regions)-1)) + 1
    ii <- round(((sections + nbands) / (2*nbands-1)) * (length(col.regions)-1)) + 1
    #ii <- sections + nbands + 1
    col <- col.regions[ii]
    for (i in seq_along(sections)) {
        section <- sections[i]
        yi <- y
        if (section < 0) {
            yi <- origin + origin - y
            section <- abs(section) - 1
        }
        baseline <- origin + section * scale
        if (all(yi <= baseline, na.rm = TRUE))
            next
        yi <- yi - baseline
        yi <- origin + pmax(pmin(yi, scale), 0)
        panel.xyarea(x, yi, border = border, col = col[i], col.line = col[i], ...)
    }
}

prepanel.horizonplot <-
    function(x, y, ..., horizonscale = NA,
             nbands = 3L,
             origin = function(y) na.omit(y)[1])
{
    if (is.function(origin))
        origin <- origin(y)
    ans <- prepanel.default.xyplot(x, y, ...)
    if (is.na(horizonscale))
        horizonscale <- max(abs(ans$ylim - origin)) / nbands
    ans$ylim <- origin + c(0, horizonscale)
    ans
}
