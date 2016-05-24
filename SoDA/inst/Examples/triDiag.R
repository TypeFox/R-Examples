triDiag <- function(diagonal, upper, lower,
                    nrow = length(diagonal), ncol = nrow) {
    value <- diag(diagonal, nrow, ncol)
    R <- row(value)
    C <- col(value)
    value[C == R + 1] <- upper
    value[C == R - 1] <- lower
    value
}
