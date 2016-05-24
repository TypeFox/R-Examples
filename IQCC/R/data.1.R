data.1 <- function(m, n, mu, Sigma)
{
    p <- dim(Sigma)[1]
    if(n == 1)
    {
        u <- matrix(nrow = m, ncol = p)
        for(i in 1:m)
        {
            N <- mvrnorm(n, mu, Sigma)
            u[i, ] <- N
        }
    }
    if(n > 1)
    {
        u <- array(dim = c(n, p, m))
        for(i in 1:m)
        {
            N <- mvrnorm(n, mu, Sigma)
            u[, , i] <- N 
        }
    }
    return(u)
}