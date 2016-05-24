## compute PCA of iris data set with svd() and sparsesvd()
library(sparsesvd)
library(Matrix)

data(iris)
M <- scale(as.matrix(iris[, 1:4]), scale=FALSE)
Ms <- Matrix(M) # not sparse, but a dMatrix

res1 <- svd(M)
res2 <- sparsesvd(Ms)

## check that eigenvalues are the same
print(res2$d, digits=3)
stopifnot(all.equal(res1$d, res2$d, tolerance=1e-12))

## these should be diagonal unit matrices
UtU <- abs(crossprod(res2$u, res1$u)) # diagonal entries may be 1 or -1
VtV <- abs(crossprod(res2$v, res1$v)) # (because sign of eigenvectors is arbitrary)
I1 <- diag(rep(1, length(res1$d)))

print(round(UtU, 12))
stopifnot(all.equal(UtU, I1, tolerance=1e-12))

print(round(VtV, 12))
stopifnot(all.equal(VtV, I1, tolerance=1e-12))

