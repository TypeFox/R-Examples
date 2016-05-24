
## not exported; for use in c.trellis and doubleYScale
mergeTrellisLegends <- function(legend, legend2, vertical = NULL)
{
    legend <- as.list(legend)
    legend2 <- as.list(legend2)
    for (space in c("top", "bottom", "left", "right")) {
        if (!is.null(legend2[[space]])) {
            if (is.null(legend[[space]])) {
                ## no conflict
                legend[[space]] <- legend2[[space]]
            } else {
                v <- vertical
                if (is.null(v))
                    v <- space %in% c("left", "right")
                legend[[space]] <-
                    list(fun = "mergedTrellisLegendGrob",
                         args = list(a = legend[[space]], b = legend2[[space]],
                         vertical = v))
            }
        }
    }
    legend <- c(legend, legend2[names(legend2) == "inside"])
    legend
}

## exported, to be called at plot time from 'legend'
mergedTrellisLegendGrob <-
    function(a, b, vertical = FALSE, border = NULL)
{
    if (is.null(a))
        return(b)
    if (is.null(b))
        return(a)
    if (!inherits(a$fun, "grob")) {
        ## fun <- a$fun
        if (is.character(a$fun)) a$fun <- as.symbol(a$fun)
        a$fun <- eval(as.call(c(a$fun, a$args)), getNamespace("lattice"))
    }
    if (!inherits(b$fun, "grob")) {
        if (is.character(b$fun)) b$fun <- as.symbol(b$fun)
        b$fun <- eval(as.call(c(b$fun, b$args)), getNamespace("lattice"))
    }
    g <- frameGrob(name = "mergedLegend")
    g <- packGrob(g, a$fun, side = if (vertical) "top" else "left", border = border)
    g <- packGrob(g, b$fun, side = if (vertical) "bottom" else "right", border = border)
    g
}

