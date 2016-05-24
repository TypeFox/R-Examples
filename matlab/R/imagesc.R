###
### $Id: imagesc.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Displays matrix C as an image with its data scaled to use the full
### color palette.
###


##-----------------------------------------------------------------------------
imagesc <- function(x = seq(ncol(C)),
                    y = seq(nrow(C)),
                    C,
                    col = jet.colors(12),
                    ...) {
    if (missing(C)) {
        if (!missing(x)) {
            if (is.null(dim(x))) {
                stop(sprintf("argument %s must be matrix-like", sQuote("x")))
            }
            C <- x
            x <- seq(ncol(C))
        } else {
            stop(sprintf("argument %s not specified", sQuote("C")))
        }
    }

    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
        stop(sprintf("increasing %s and %s values expected",
                     sQuote("x"), sQuote("y")))
    }

    if (!(is.numeric(C) && is.matrix(C))) {
        stop(sprintf("argument %s must be matrix", sQuote("C")))
    }

    range.C <- range(C)
    min.lim.C <- round(range.C[1])
    max.lim.C <- round(range.C[2])

    graphics::image(x    = x,
                    y    = y,
                    z    = t(C)[,nrow(C):1],
                    zlim = c(min.lim.C, max.lim.C),
                    axes = FALSE,
                    col  = col,
                    ...)

    pretty.axp <- function(axp.name = c("xaxp", "yaxp"), ...) {
        axp.name <- match.arg(axp.name)

        dots <- list(...)
        axp <- if (axp.name %in% names(dots)) {
                  dots[[axp.name]]
               } else {
                  par(axp.name)
               }

        pretty(axp[1]:axp[2], axp[3])
    }

    at.x <- pretty.axp("xaxp", ...)
    axis(SIDE.BELOW, at.x)
    at.y <- pretty.axp("yaxp", ...)
    axis(SIDE.LEFT, at.y, labels = rev(at.y), las = LAS.HORIZONTAL)
    box()

    invisible(NULL)
}

