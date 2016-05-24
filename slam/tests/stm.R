
library("slam")
set.seed(20090626)

###

x <- sample(0:5, 100, T, prob=c(.8,rep(.04,5)))
x <- matrix(as.logical(x), nrow = 20,
     dimnames = list(rows = 1:20, cols = LETTERS[1:5]))
x

xst <- as.simple_triplet_matrix(x)
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

xnast <- as.simple_triplet_matrix(xna)
xnast

## default method
identical(rowSums(xna), row_sums(xna))
identical(colSums(xna), col_sums(xna))
identical(rowMeans(xna), row_means(xna))
identical(colMeans(xna), col_means(xna))

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
identical(tcrossprod(x[1:10,], x[11:20,]),
	  tcrossprod_simple_triplet_matrix(xst[1:10,], xst[11:20,]))

x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 3)
x
s <- as.simple_triplet_matrix(x)

identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(s))
identical(tcrossprod(x), tcrossprod_simple_triplet_matrix(s, x))
identical(tcrossprod(x[2:3,], x[1,, drop = FALSE]), 
	  tcrossprod_simple_triplet_matrix(s[2:3,], s[1,]))
identical(tcrossprod(x[1,, drop = FALSE], x[2:3,]), 
	  tcrossprod_simple_triplet_matrix(s[1,], s[2:3,]))

###


