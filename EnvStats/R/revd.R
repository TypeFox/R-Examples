revd <-
function (n, location = 0, scale = 1) 
{
    na.index <- is.na(scale)
    if (any(!na.index)) 
        if (scale[!na.index] < .Machine$double.eps) 
            stop("All values of 'scale' must be positive.")
    location - scale * log(-log(runif(n)))
}
