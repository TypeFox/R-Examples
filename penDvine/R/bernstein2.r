bernstein2 <- function(x,v,n)
{
    return((choose(n, v) * x^v * (1 - x)^(n - v)) * (n + 1))
}
