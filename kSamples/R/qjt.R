qjt <- function (p, nn) 
{
    dist <- Harding(nn)
    if (is.nan(dist[1])) 
        stop("can't compute due to numerical instability\n")
    nn <- sort(as.integer(nn))
    nvec <- rev(cumsum(rev(nn)))
    k <- length(nn)
    L <- sum(nn[1:(k - 1)] * nvec[2:k])
    dist <- c(0, cumsum(dist))
    x <- c(-Inf, 0:L)
    out <- findInterval(p, dist)
    out[p %in% dist] <- out[p %in% dist] - 1
    x[out + 1]
}
