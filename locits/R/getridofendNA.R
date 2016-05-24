getridofendNA <-
function (x) 
{
    x[is.na(x)] <- 0
    return(x)
    n <- length(x)
    halfn <- n/2
    x1 <- x[1:halfn]
    x2 <- x[(halfn + 1):n]
    x2 <- rev(x2)
    ix1 <- which(is.na(x1))
    if (length(ix1 > 0)) 
        x1[ix1] <- rev(x1[(max(ix1) + 1):(max(ix1) + length(ix1))])
    ix2 <- which(is.na(x2))
    if (length(ix2 > 0)) 
        x2[ix2] <- rev(x2[(max(ix2) + 1):(max(ix2) + length(ix2))])
    c(x1, rev(x2))
}
