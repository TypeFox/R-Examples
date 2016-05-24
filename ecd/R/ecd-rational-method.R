#' Utility to convert a numeric to a rational
#' 
#' Convert a numeric x to rational p/q, which is then used for polynomial construction
#'
#' @param x numeric
#'
#' @return a vector of two integers, representing numerator and denominator
#'
#' @keywords solve
#'
#' @export 
#'
#' @examples
#' pq <- ecd.rational(2.5)
### <======================================================================>
ecd.rational <- function(x)
{
    if (length(x) != 1 | !is.finite(x)) {
        stop("Input must be length-one finite numeric!")
    }
    
    cycles = 10
    max.denominator = 100
    
    a0 <- rep(0, length(x))
    b0 <- rep(1, length(x))
    A <- matrix(b0)
    B <- matrix(floor(x))
    r <- as.vector(x) - drop(B)
    
    len <- 0
    while(any(r > 1/max.denominator) && (len <- len + 1) <= cycles) {
        a <- a0
        b <- b0
        a <- 1
        r <- 1/r
        b <- floor(r)
        r <- r - b
        A <- cbind(A, a)
        B <- cbind(B, b)
    }
    
    pq1 <- cbind(b0, a0)
    pq <- cbind(B[, 1], b0)
    len <- 1
    while((len <- len + 1) <= ncol(B)) {
        pq0 <- pq1
        pq1 <- pq
        pq <- B[, len] * pq1 + A[, len] * pq0
    }
    unname(as.vector(pq))
}