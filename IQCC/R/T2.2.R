T2.2 <- function(datum2, estat, n)
{
    if(n == 1)
    {
        T2 <- (t(datum2 - estat[[1]]) %*% solve(estat[[2]]) %*% (datum2 - estat[[1]]))
    }
    if(n > 1)
    {
        media <- colMeans(datum2)
        T2 <- n * (t(media - estat[[1]]) %*% solve(estat[[2]]) %*% (media - estat[[1]]))
    }
    return(T2)
}