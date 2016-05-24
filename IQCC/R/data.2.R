data.2 <- function(estat, n, delta = 0, p)
{
    N <- mvrnorm(n, c(estat[[1]]) + delta, matrix(c(estat[[2]]), p, p))
}