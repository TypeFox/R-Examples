triDiag3 <- function(diagonal, upper, lower,
                    nrow = length(diagonal), ncol = nrow) {
    value <- matrix(0, nrow = nrow, ncol = ncol)
    r <-max(nrow, ncol)
    if(r > 1) {
        nu <- min(nrow, ncol-1)
        nl <- min(nrow-1, ncol)
        index <- outer((0:nu)*(nrow+1), 0:2, `+`)
        value[index[1:min(nrow, ncol), 2]] <- diagonal
        if(nu > 0)
            value[index[-1, 1]] <- upper
        if(nl > 0)
            value[index[1:nl, 3]] <- lower
    }
    value
}
