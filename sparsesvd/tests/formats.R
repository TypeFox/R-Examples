## any sparse matrix format that inherits from dMatrix should work
library(sparsesvd)
library(Matrix)

M <- rbind(
  c(20, 10, 15,  0,  2),
  c(10,  5,  8,  1,  0),
  c( 0,  1,  2,  6,  3),
  c( 1,  0,  0, 10,  5))

res1 <- sparsesvd(as(M, "dgCMatrix")) # standard format (column-compressed)
res2 <- sparsesvd(as(M, "dgeMatrix")) # dense matrix
res3 <- sparsesvd(as(M, "dgTMatrix")) # triple format
## -- row-compressed form cannot be converted to dgCMatrix
## res4 <- sparsesvd(as(M, "dgRMatrix")) # row-compressed

## check that eigenvalues are the same
stopifnot(all.equal(res1$d, res2$d, tolerance=1e-12))
stopifnot(all.equal(res1$d, res3$d, tolerance=1e-12))
## stopifnot(all.equal(res1$d, res4$d, tolerance=1e-12))

## special classes for symmetric matrices
A <- crossprod(M)
res1a <- sparsesvd(as(A, "dgCMatrix")) # standard format (column-compressed)
## -- symmetric matrix can only be converted if already in row-compressed form
## res2a <- sparsesvd(as(A, "dsTMatrix")) # symmetric triplet format
res3a <- sparsesvd(as(A, "dsCMatrix")) # symmetric column-compressed

## check that eigenvalues are the same and consistent with M
stopifnot(all.equal(res1a$d, (res1$d)^2, tolerance=1e-12))
## stopifnot(all.equal(res1a$d, res2a$d, tolerance=1e-12))
stopifnot(all.equal(res1a$d, res3a$d, tolerance=1e-12))
