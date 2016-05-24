rhypoexp <- function(n = 1L, rate = 1.0) {

    if ((rlen <- length(rate)) == 1L) {
        return(rexp(n = n, rate = rate))   
    }

    if ((nlen <- length(n)) > 1L) {
        n <- nlen
    }

#     rowSums(matrix(vapply(X         = rate,
#                           FUN       = rexp,
#                           FUN.VALUE = numeric(n),
#                           n         = n),
#                    nrow = n,
#                    ncol = length(rate)))

    colSums(matrix(data = rexp(n    = n * rlen,
                               rate = rate),
                   nrow = rlen,
                   ncol = n))
}
