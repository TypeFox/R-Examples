LocalCoarsen <- function(x, w, IsKnot){
    n <- length(x)
    K <- (1:n) * IsKnot
    K <- K[K > 0]
    x2 <- x[K]
    w2 <- w[K]
    for (k in 1:(length(K) - 1)){
        if (K[k + 1] > (K[k] + 1)){
            ind <- (K[k] + 1):(K[k + 1] - 1)
            lambda <- (x[ind] - x2[k])/(x2[k + 1] - x2[k])
            w2[k] <- w2[k] + sum(w[ind] * (1 - lambda))
            w2[k + 1] <- w2[k + 1] + sum(w[ind] * lambda)
        }
    }
    w2 <- w2 / sum(w2)
    return(list(x2 = matrix(x2, ncol = 1), w2 = matrix(w2, ncol = 1)))
}
