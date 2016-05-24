LocalExtend <- function(x, IsKnot, x2, phi2){
    n <- length(x)
    K <- (1:n) * IsKnot
    K <- K[K > 0]
    phi <- 1:n * 0
    phi[K] <- phi2
    for (k in 1:(length(K) - 1)){
        if (K[k + 1] > (K[k] + 1)){
            ind <- (K[k] + 1):(K[k + 1] - 1)
            lambda <- (x[ind] - x2[k])/(x2[k + 1] - x2[k])
            phi[ind] <- (1 - lambda) * phi2[k] + lambda * phi2[k + 1]
        }
    }
    return(matrix(phi, ncol = 1))
}
