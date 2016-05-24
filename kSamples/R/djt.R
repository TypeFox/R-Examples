djt <- function (x, nn) 
{
    dist <- Harding(nn)
    if (is.nan(dist[1])) 
        stop("can't compute due to numerical instability\n")
    nn <- sort(as.integer(nn))
    nvec <- rev(cumsum(rev(nn)))
    k <- length(nn)
    L <- sum(nn[1:(k - 1)] * nvec[2:k])
    pos <- match(x, 0:L)
    d <- numeric(length(x))
    d[!is.na(pos)] <- dist[pos[!is.na(pos)]]
    d
}
