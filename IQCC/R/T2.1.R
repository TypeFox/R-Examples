T2.1 <- function(estat, m, n)
{
    t2 <- vector()
    if(n == 1)
    {
        for (i in 1:m)
        {
            T2 <- (t(estat[[3]][i, ]) %*% solve(estat[[2]]) %*% (estat[[3]][i, ]))
            t2 <- c(t2, T2)
        }
    }
    if(n > 1)
    {
        for (i in 1:m)
        {
            T2 <- n * (t(estat[[3]][i, ] - estat[[1]]) %*% solve(estat[[2]]) %*% (estat[[3]][i, ] - estat[[1]]))
            t2 <- c(t2, T2)
        }
    }
    return(t2)
}