ANDP <- function(b, bE) {
    if (!is.vector(b) | !is.vector(bE)) 
        stop("Function is only defined for vectors")
    if (length(b) != length(bE)) 
        stop("Function is only defined for vectors of the same length")
    if (length(b) < 1) 
        stop("Function is only defined for non-empty vectors")
    
    n <- length(b)
    Nf <- sum((abs(b[2:n] - b[1:(n - 1)]) != 0)) + 1
    NfE <- sum((abs(bE[2:n] - bE[1:(n - 1)]) != 0)) + 1
    s <- abs(Nf - NfE)
    return(s)
} 
