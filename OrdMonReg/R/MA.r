MA <- function(g, w, A = NA, a, b){

    ## Compute the value of M(A) as defined in Barlow, Bartholomew, 
    ## Bremner & Brunk (1972), p. 57
    if (max(is.na(A)) == TRUE){A <- 1:length(g)}
    
    a.up <- max(a[A], na.rm = TRUE)
    b.low <- min(b[A], na.rm = TRUE)
    
    res <- NA
    if (a.up <= b.low){res <- max(min(sum((g * w)[A]) / sum(w[A]), b.low), a.up)}
    return(res)
}
