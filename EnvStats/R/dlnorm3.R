dlnorm3 <-
function (x, meanlog = 0, sdlog = 1, threshold = 0) 
{
    y <- dlnorm(x = x - threshold, meanlog = meanlog, sdlog = sdlog)
    if (!is.null(Names <- names(x))) 
        names(y) <- rep(Names, length = length(y))
    y
}
