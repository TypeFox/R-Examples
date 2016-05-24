PCE <- function(b, bE) {
    if (!is.vector(b) | !is.vector(bE)) 
        stop("Function is only defined for vectors")
    if (length(b) != length(bE)) 
        stop("Function is only defined for vectors of the same length")
    if (length(b) < 1) 
        stop("Function is only defined for non-empty vectors")
    
    n <- length(b)
    s <- 100 * sum(b != bE)/n
    return(s)
} 
