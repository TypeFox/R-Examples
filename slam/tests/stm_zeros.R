
library("slam")
set.seed(20091012)


###
as.simple_triplet_matrix_zeros <-
function(x) {
    x <- list(
        i        = rep(seq_len(nrow(x)), ncol(x)),
        j        = rep(seq_len(ncol(x)), each = nrow(x)),
        v        = c(x),
        nrow     = nrow(x),
        ncol     = ncol(x),
        dimnames = dimnames(x)
    )
    class(x) <-  "simple_triplet_matrix"
    x
}

x <- sample(0:5, 100, T, prob=c(.8,rep(.04,5)))
x <- matrix(as.logical(x), nrow = 20,
     dimnames = list(rows = 1:20, cols = LETTERS[1:5]))
x

xst <- as.simple_triplet_matrix_zeros(x)
xst

identical(rowSums(x), row_sums(xst))
identical(colSums(x), col_sums(xst))
identical(rowMeans(x), row_means(xst))
identical(colMeans(x), col_means(xst))

## NAs

xna <- x
n <- prod(dim(x))
is.na(xna) <- sample(seq_len(n), ceiling(n * .1))
xna

xnast <- as.simple_triplet_matrix_zeros(xna)
xnast

identical(rowSums(xna), row_sums(xnast))
identical(colSums(xna), col_sums(xnast))
identical(rowMeans(xna), row_means(xnast))
identical(colMeans(xna), col_means(xnast))

identical(rowSums(xna, na.rm = TRUE), row_sums(xnast, na.rm = TRUE))
identical(colSums(xna, na.rm = TRUE), col_sums(xnast, na.rm = TRUE))
identical(rowMeans(xna, na.rm = TRUE), row_means(xnast, na.rm = TRUE))
identical(colMeans(xna, na.rm = TRUE), col_means(xnast, na.rm = TRUE))

## cross-product

identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(xst))
identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(xst, x))

x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 3)
x
s <- as.simple_triplet_matrix_zeros(x)

identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(s))
identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(s, x))

##
identical(as.matrix(s * x), x * x)
identical(as.matrix(x * s), x * x)
identical(as.matrix(s * s), x * x)

identical(as.matrix(s + s), x + x)

###


