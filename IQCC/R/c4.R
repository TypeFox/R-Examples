c4 <- function(n)
{
    c <- (sqrt(2 / (n-1))) * (gamma(n/2) / gamma((n-1) / 2))
    return(c)
}