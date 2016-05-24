## extra triDiag functions

triDiag2S <- function(diagonal, upper, lower, r = length(diagonal)) {
    value <- diag(diagonal, r)
    if(r > 1) {
        rseq <- 1:(r-1)
        index <- cbind(rseq+1, rseq)
        value[index] <- lower
        index <- cbind(rseq, rseq+1)
        value[index] <- upper
    }
    value
}

triDiag3S <- function(diagonal, upper, lower,
                      r = length(diagonal)) {
    value <- matrix(0, r, r)
    if(r > 0) {
        r1 <- r-1
        index <- outer((0:r1)*(r+1), 0:2, "+")
        value[index[,2]] <- diagonal
        value[index[-1, 1]] <- upper
        value[index[-r, 3]] <- lower
    }
    value
}
