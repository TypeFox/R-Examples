MNChPD <- function(b, bE) {
    if (!is.vector(b) | !is.vector(bE)) 
        stop("Function is only defined for vectors")
    if (length(b) != length(bE)) 
        stop("Function is only defined for vectors of the same length")
    if (length(b) < 1) 
        stop("Function is only defined for non-empty vectors")
    n <- length(b)
    
    CPs_b <- c(0, 2^b[2:n] - 2^b[1:(n - 1)])
    CPs <- c(0, 2^bE[2:n] - 2^bE[1:(n - 1)])
    ps_b <- which(CPs_b != 0)
    
    if (length(ps_b) != 0) {
        sd <- 0
        for (i in 1:length(ps_b)) {
            k <- ps_b[i]
            tipo = 2^b[k] - 2^b[k - 1]
            ps <- which(CPs == tipo)
            d <- min(abs(ps - k))
            sd <- sd + d
        }
        s <- sd/length(ps_b)
    } else {
        s <- 0
    }
    return(s)
} 
