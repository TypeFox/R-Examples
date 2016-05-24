"image.asc" <- function (x, col = gray((240:1)/256), clfac = NULL, ...)
{
    ## Verifications
    if (!inherits(x, "asc"))
        stop("not an \"asc\" object")

    ## Coordinates of the pixels
    xy <- getXYcoords(x)
    xx <- xy$x
    yy <- xy$y

    ## If the variable is numeric
    if (attr(x, "type") == "numeric")
        image(x = xx, y = yy, x, asp = 1, col = col, ...)

    ## For a factor: creates colors
    if (attr(x, "type") == "factor") {
        if (is.null(clfac)) {
            clfac <- rainbow(nlevels(x))
            clfac <- clfac[as.numeric(levels(factor(x)))]
        }
        image(x = xx, y = yy, x, asp = 1, col = clfac, ...)
    }
}
