##
##  p ol y p o w . R  Polynomial Powers
##


polypow <- function(p, n){
    if ( !is.vector(p, mode="numeric") && !is.vector(p, mode="complex") )
        stop("Arguments 'p' must be a real or complex vector.")
    if ( !is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 0 )
        stop("Argument 'n' must be a non-negative integer.")

    pp <- c(1)
    while (n > 0) {
    	pp <- polymul(pp, p)
    	n <- n - 1
    }

    return(pp)
}


polytrans <- function(p, q){
    if ( (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex")) ||
         (!is.vector(q, mode="numeric") && !is.vector(q, mode="complex")) )
        stop("Arguments 'p' and 'q' must be real or complex vectors.")

    n <- length(p)
    if (length(p) == 1)
        return(p)

    pt <- 0
    for (i in 1:n) {
    	pt <- polyadd(pt, p[i]*polypow(q, n-i))
    }

    return(pt)
}
