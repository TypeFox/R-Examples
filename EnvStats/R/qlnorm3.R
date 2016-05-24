qlnorm3 <-
function (p, meanlog = 0, sdlog = 1, threshold = 0) 
{
    q <- qlnorm(p = p, meanlog = meanlog, sdlog = sdlog)
    lq <- length(q)
    q <- q + (threshold <- rep(threshold, len = lq))
    if (any(index <- rep(p, length = lq) == 0 & !is.na(q))) 
        q[index] <- threshold[index]
    if (!is.null(Names <- names(p))) 
        names(q) <- rep(Names, length = lq)
    q
}
