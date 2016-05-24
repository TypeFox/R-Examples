pevd <-
function (q, location = 0, scale = 1) 
{
    na.index <- is.na(scale)
    if (any(!na.index)) 
        if (any(scale[!na.index] < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
    p <- exp(-exp(-(q - location)/scale))
    if (!is.null(Names <- names(q))) 
        names(p) <- rep(Names, length = length(p))
    p
}
