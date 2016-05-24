
## CB 2009/5,6,10 2010/6 2013/10


## NOTE the C code does not use long double for accumulation.
.means_simple_triplet_matrix <-
function(x, DIM, na.rm)
{
    s <- .Call(R_sums_stm, x, DIM, na.rm)
    n <- c(x$nrow, x$ncol)[-DIM]
    if (na.rm) {
	x$v <- is.na(x$v)
	nna <- .Call(R_sums_stm, x, DIM, FALSE)
	s / (n - nna)
    }
    else
	s /  n
}


## R interfaces

row_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("row_sums")

row_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::rowSums(x, na.rm, dims, ...)

row_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .Call(R_sums_stm, x, 1L, na.rm)

row_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)
row_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)

col_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("col_sums")

col_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::colSums(x, na.rm, dims, ...)

col_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .Call(R_sums_stm, x, 2L, na.rm)

col_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)
col_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)

row_means <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("row_means")

row_means.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::rowMeans(x, na.rm, dims, ...)

row_means.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 1L, na.rm)

row_means.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowMeans(x, na.rm = na.rm, dims = dims, ...)
row_means.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowMeans(x, na.rm = na.rm, dims = dims, ...)

col_means <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("col_means")

col_means.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::colMeans(x, na.rm, dims, ...)

col_means.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    .means_simple_triplet_matrix(x, DIM = 2L, na.rm)

col_means.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colMeans(x, na.rm = na.rm, dims = dims, ...)
col_means.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colMeans(x, na.rm = na.rm, dims = dims, ...)

row_norms <-
function(x, p = 2)
{
    if(p == 2)
        sqrt(row_sums(x ^ 2))
    else if(p == 1)
        row_sums(abs(x))
    else if(p == Inf)
        c(rollup(abs(x), 2L, FUN = max))
    else
        row_sums(abs(x) ^ p) ^ (1/p)
}

col_norms <-
function(x, p = 2)
{
    if(p == 2)
        sqrt(col_sums(x ^ 2))
    else if(p == 1)
        col_sums(abs(x))
    else if(p == Inf)
        c(rollup(abs(x), 1L, FUN = max))
    else
        col_sums(abs(x) ^ p) ^ (1/p)
}

##
.nnzero <- 
function(x, scale = FALSE) {
    v <- c("simple_triplet_matrix", "simple_sparse_array")
    if (inherits(x, v))
	v <- x$v
    else {
	x <- as.array(x)
	v <- x
    }
    v <- v == vector(typeof(v), 1L)
    v <- v + 1L
    n <- length(v)
    v <- tabulate(v, 2L)
    v <- c(v, n - sum(v))
    names(v) <- c("nnzero", "nzero", NA)
    if (scale)
	v <- v / prod(dim(x))
    v
}

###
