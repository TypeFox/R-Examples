## Simple functions for converting "matrix" type objects into the
## sparse "column major order" (CSC, modulo offsets) format used by
## SYMPHONY.

## matind: vector of the row indices corresponding to each entry of
##   value 
## values: vector of the values of nonzero entries of the constraint
##   matrix in column order.

make_csc_matrix <-
function(x)
    UseMethod("make_csc_matrix")

make_csc_matrix.matrix <-
function(x)
{
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")
   
    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1] - 1L,
         values = x[ind])
}

make_csc_matrix.simple_triplet_matrix <-
function(x)
{
    if(!inherits(x, "simple_triplet_matrix"))
        stop("Argument 'x' must be of class 'simple_triplet_matrix'.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}
