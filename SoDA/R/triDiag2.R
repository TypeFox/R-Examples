triDiag2 <- function(diagonal, upper, lower,
                    nrow = length(diagonal), ncol = nrow) {
    value <- diag(diagonal, nrow = nrow, ncol = ncol)
    n <- min(nrow, ncol-1)
    if(n>0) {
        rseq <- 1:n
        value[cbind(rseq, rseq+1)] <- upper
    }
    n <- min(nrow-1, ncol)
    if(n > 0) {
        rseq <- 1:n
        value[cbind(rseq+1, rseq)] <- lower
    }
    value
}
