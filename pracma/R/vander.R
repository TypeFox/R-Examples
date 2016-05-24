###
### VANDER.R  Vandermonde matrix
###


vander <- function(x) {
    n <- length(x)
    if (n == 0) return(matrix(0, nrow=0, ncol=0))
    if ((!is.numeric(x) && !is.complex(x)) || is.array(x))
        stop("Argument 'x' must be a numeric or complex  vector.")

    A <- outer(x, seq(n-1, 0), "^")
    return(A)
}
