minK1 <- function(g, w, n){
    r1 <- BoundedIsoMean(g[1:n], w[1:n])
    r2 <- g[(n + 1):(2 * n)]
    res <- c(r1, r2)
    return(res)
    }
