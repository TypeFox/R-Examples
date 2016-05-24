##
##  b e r n o u l l i . R  Bernoulli numbers and polynomials
##


bernoulli <- function(n, x) {
    stopifnot(length(n) == 1, floor(n) == ceiling(n), n >= 0)
    if (missing(x)) {
        if (n == 0) return(1.0)
        bf <- function(n) -n * zeta(1-n)
        bx <- c(1.0, sapply(1:n, bf))
        bx[2] <- -0.5
    } else {
        stopifnot(is.numeric(x))
        bn <- bernoulli(n)
        cn <- choose(n, 0:n)
        bx <-polyval(cn*bn, x)
    }
    bx
}


# euler <- function(n) {
#     stopifnot(length(n) >= 1, floor(n) == ceiling(n))
#     if (n == 0) return(1)
#     en <- numeric(n)
#     en[1] <- 0.0
#     en[2] <- -1.0
#     if (n < 4) return(c(1.0, en[1:n]))
#     for (m in 2:floor(n/2)) {
#         s <- 1.0
#         for (k in 1:(m-1)) {
#             r = 1.0
#             for (j in 1:(2*k)) {
#                 r <- r * (2*m - 2*k + j) / j
#             }
#             s <- s + r*en[2*k]
#         }
#         en[2*m] <- (-s)
#     }
#     c(1.0, en)
# }
