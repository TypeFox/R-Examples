dPval <- function(p, u, theta, sigma2){
    
    eta <- sqrt(u ^ 2 + sigma2)
    z <- qnorm(p / 2)
       
    if (1 == 0){
    ## using folded normal
    t1 <- dnorm(-u * z, mean = theta, sd = eta) + dnorm(u * z, mean = theta, sd = eta)
    t2 <- dnorm(-u * z, mean = 0, sd = u) + dnorm(u * z, mean = 0, sd = u)
    }
    
    ## simplified
    t1 <- dnorm((-u * z - theta) / eta) / eta + dnorm((u * z - theta) / eta) / eta
    t2 <- 2 * dnorm(z) / u
    
    res <- t1 / t2
    return(res)
}

