devd <-
function (x, location = 0, scale = 1) 
{
    na.index <- is.na(scale)
    if (any(!na.index)) 
        if (any(scale[!na.index] < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
    z <- (x - location)/scale
    y <- (1/scale) * exp(-z - exp(-z))
    if (!is.null(Names <- names(x))) 
        names(y) <- rep(Names, length = length(y))
    y
}
