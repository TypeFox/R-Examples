TD <- function(b, bE, t) {
    if (!is.vector(b) | !is.vector(bE)) 
        stop("Function is only defined for vectors")
    if (length(b) != length(bE)) 
        stop("Function is only defined for vectors of the same length")
    if (length(b) < 1) 
        stop("Function is only defined for non-empty vectors")
    if (t <= 0) 
        stop("Tolerance value 't' should be positive")
    
    n <- length(b)
    CPs_b <- c(0, 2^b[2:n] - 2^b[1:(n - 1)])
    CPs <- c(0, 2^bE[2:n] - 2^bE[1:(n - 1)])
    ps_b <- which(CPs_b != 0)
    if (length(ps_b) != 0) {
        sd <- numeric()
        for (i in 1:length(ps_b)) {
            k <- ps_b[i]
            tipo = 2^b[k] - 2^b[k - 1]
            ps <- which(sign(CPs) == sign(tipo))
            d <- min(abs(ps - k))
            if (d < t) 
                sd <- c(sd, d)
        }
    if(length(sd)>1){
        s <- sqrt(sum(sd^2)/length(sd))
    } else {
       s <- NA
       warning("Can't compute TD, tolerance value too low.")
    }
    } else {
        s <- 0
    }
    return(s)
} 
