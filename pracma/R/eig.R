###
### EIG.R  Eigenvalues
###


eig <- function(a) {
    if (length(a) == 0) return(matrix(0, nrow=0, ncol=0))
    if (length(a) == 1) return(a)
    if ((!is.numeric(a) && !is.complex(a)) || !is.matrix(a))
        stop("Argument 'a' must be a numeric or complex matrix.")
    if (nrow(a) != ncol(a))
        stop("Matrix 'a' must be square matrix.")

    eigen(a, only.values=TRUE)$values
}
