minK2 <- function(g, w, n){
    r1 <- g[1:n]
    i <- (n + 1):(2 * n)
    r2 <- BoundedIsoMean(g[i], w[i])
    res <- c(r1, r2)
    return(res)
    }
