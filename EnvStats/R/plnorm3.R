plnorm3 <-
function (q, meanlog = 0, sdlog = 1, threshold = 0) 
{
    p <- plnorm(q = q - threshold, meanlog = meanlog, sdlog = sdlog)
    if (!is.null(Names <- names(q))) 
        names(p) <- rep(Names, length = length(p))
    p
}
