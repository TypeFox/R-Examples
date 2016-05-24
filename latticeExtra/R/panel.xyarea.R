##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

panel.xyarea <- function(x, ...)
    UseMethod("panel.xyarea")

## Plot a series as a filled polygon connected at given origin (on y axis).
## With groups, acts like panel.superpose, but with polygon style settings.
panel.xyarea.default <-
    function(x, y, groups = NULL, origin = NULL, horizontal = FALSE,
             col = if (is.null(groups)) plot.polygon$col else superpose.polygon$col,
             col.line = if (is.null(groups)) plot.polygon$col else superpose.polygon$col,
             border = if (is.null(groups)) plot.polygon$border else superpose.polygon$border,
             lty = if (is.null(groups)) plot.polygon$lty else superpose.polygon$lty,
             lwd = if (is.null(groups)) plot.polygon$lwd else superpose.polygon$lwd,
             alpha = if (is.null(groups)) plot.polygon$alpha else superpose.polygon$alpha,
             ..., fill, panel.groups = panel.xyarea)
{
    plot.polygon <- trellis.par.get("plot.polygon")
    superpose.polygon <- trellis.par.get("superpose.polygon")
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) == 0) return()
    if (!is.null(groups)) {
        ## NOTE superpose does not handle 'border' argument, so pass it as col.line
        panel.superpose(x, y, ..., groups = groups, panel.groups = panel.groups,
                        col = col, col.line = col.line, lty = lty, lwd = lwd, border = border,
                        alpha = alpha, origin = origin, horizontal = horizontal)
    } else {
        if (!missing(col.line))
            col <- col.line
        if (horizontal == TRUE) {
            ## actually means origin is vertical. for consistency with panel.xyplot.
            xlim <- current.panel.limits()$xlim
            if (is.null(origin))
                origin <- xlim[1]
            infi <- is.infinite(x)
            x[infi] <- ifelse(x[infi] > 0, max(xlim), min(xlim))
        } else {
            ## default case; origin is horizontal
            ylim <- current.panel.limits()$ylim
            if (is.null(origin))
                origin <- ylim[1]
            infi <- is.infinite(y)
            y[infi] <- ifelse(y[infi] > 0, max(ylim), min(ylim))
        }
        stopifnot(is.numeric(origin))
        ## need to split up the series into chunks without any missing values
        ## (because NAs split the polygon)
        xy <- data.frame(x = x, y = y)
        ## order by ordinate values
        ord <- if (horizontal) order(xy$y) else order(xy$x)
        xy <- xy[ord,]
        ok <- complete.cases(xy)
        runs <- rle(ok)
        ## assign unique values to each chunk, and NAs between (dropped by 'split')
        runs$values[runs$values == TRUE] <- seq_len(sum(runs$values))
        runs$values[runs$values == FALSE] <- NA
        ## expand into long format
        chunks <- inverse.rle(runs)
        lapply(split(xy, chunks), function(xyi, ...) {
            x <- xyi$x
            y <- xyi$y
            ## drop ends of series to the origin; the polygon will be joined up at that level
            if (horizontal == TRUE) {
                ## non-default case
                yy <- c(head(y,1), y, tail(y,1))
                xx <- c(origin, x, origin)
            } else {
                ## default case
                xx <- c(head(x,1), x, tail(x,1))
                yy <- c(origin, y, origin)
            }
            ## we need to catch the 'fill' argument from panel.superpose, otherwise over-rides 'col'
            panel.polygon(xx, yy, alpha = alpha, col = col, border = border, lty = lty, lwd = lwd, ...)
        }, ...)
    }
}

panel.xyarea.ts <- function(x, y = x, ...)
{
    panel.xyarea(as.vector(time(x)), y, ...)
}

panel.xyarea.zoo <- function(x, y = x, ...)
{
    panel.xyarea(zoo::index(x), zoo::coredata(y), ...)
}

## A slightly modified copy of panel.qqmath
panel.qqmath.xyarea <-
    function(x, y = NULL,
             f.value = NULL,
             distribution = qnorm,
             qtype = 7,
             groups = NULL, ...,
             tails.n = 0)
{
    x <- as.numeric(x)
    distribution <-
        if (is.function(distribution)) distribution 
        else if (is.character(distribution)) get(distribution)
        else eval(distribution)
    nobs <- sum(!is.na(x))
    if (!is.null(groups))
        panel.xyarea(x, y = NULL,
                     f.value = f.value,
                     distribution = distribution,
                     qtype = qtype,
                     groups = groups,
                     panel.groups = panel.qqmath.xyarea,
                     ...,
                     tails.n = tails.n)
    else if (nobs)
    {
        if (is.null(f.value)) # exact data instead of quantiles
        {
            panel.xyarea(x = distribution(ppoints(nobs)),
                         y = sort(x),
                         ...)
        }
        else
        {
            pp <- if (is.numeric(f.value)) f.value else f.value(nobs)
            if (tails.n > 0)
            {
                ## use exact data for tails of distribution
                tails.n <- min(tails.n, nobs %/% 2)
                ppd <- ppoints(nobs)
                ## omit probabilities within the exact tails
                pp <- pp[(pp > ppd[tails.n] &
                          pp < ppd[nobs + 1 - tails.n])]
                ## add on probs corresponding to exact tails
                pp <- c(head(ppd, tails.n), pp, tail(ppd, tails.n))
                ## must use a quantile type that recovers exact values:
                qtype <- 1
            }
            xx <- distribution(pp)
            yy <- quantile(x, pp, 
                                  names = FALSE,
                                  type = qtype,
                                  na.rm = TRUE)
            panel.xyarea(x = xx, y = yy, ...)
        }
    }
}
