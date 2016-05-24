##
##  p ol y m u l .R  Polynomial Multiplication
##


polymul <- function(p, q){
    if ( (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex")) ||
         (!is.vector(q, mode="numeric") && !is.vector(q, mode="complex")) )
        stop("Arguments 'p' and 'q' must be real or complex vectors.")

    n <- length(p); m <- length(q)
    if (n <= 1 || m <= 1) return(p*q)

    r <- rep(0, n+m-1)
    for (i in seq(along=q)) {
        r <- r + c(rep(0, i-1), p * q[i], rep(0, m-i))
    }

    while (r[1] == 0 && length(r) > 1)
        r <- r[2:length(r)]

    return(r)
}
