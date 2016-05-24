Hij <- function(theta, sigma, y, u, teststat){
       
    n <- length(y)
    k <- 1 + floor(n / 2)
    eta <- sqrt(u ^ 2 + sigma ^ 2)

    ## compute Hij
    Hij <- matrix(NA, nrow = n, ncol = k)
    i <- 1:n

    Hij[, 1] <- pnorm((u[i] * teststat[2] - theta) / eta[i]) - pnorm((-u[i] * teststat[2] - theta) / eta[i])
    for (j in 2:(k - 1)){
        Hij[, j] <- pnorm(( teststat[2 * j    ] * u[i] - theta) / eta[i]) -
                    pnorm(( teststat[2 * j - 2] * u[i] - theta) / eta[i]) +
                    pnorm((-teststat[2 * j - 2] * u[i] - theta) / eta[i]) -
                    pnorm((-teststat[2 * j    ] * u[i] - theta) / eta[i])
    }
    Hij[, k] <- 1 - pnorm((teststat[2 * k - 2] * u[i] - theta) / eta[i]) + pnorm((- teststat[2 * k - 2] * u[i] - theta) / eta[i])
        
    return(list("Hij" = Hij))    
}
