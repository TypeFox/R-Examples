`orglpoints` <-
    function (object, display = "sites", choices = 1:3, radius, col = "black",
              ...)
{
    x <- scores(object, display = display, choices = choices, ...)
    ## default radius
    if (missing(radius))
        radius <- max(apply(x, 2, function(z) diff(range(z))))/100
    ## honor cex
    cex <- match.call(expand.dots = FALSE)$...$cex
    if (!is.null(cex))
        radius <- cex * radius
    ## make a color vector, handle factors
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nrow(x))
    rgl.spheres(x, radius = radius, col = col, ...)
    invisible()
}

