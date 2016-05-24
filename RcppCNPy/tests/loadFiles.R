#!/usr/bin/Rscript

library(RcppCNPy)

fmat <- npyLoad("fmat.npy")
print(fmat)
stopifnot(all.equal(fmat, t(matrix(seq(0,11) * 1.1, 4, 3))))

fmat <- npyLoad("fmat.npy.gz")
print(fmat)
stopifnot(all.equal(fmat, t(matrix(seq(0,11) * 1.1, 4, 3))))

imat <- npyLoad("imat.npy", "integer")
print(imat)
stopifnot(all.equal(imat, t(matrix(seq(0,11), 4, 3))))

fvec <- npyLoad("fvec.npy")
print(fvec)
stopifnot(all.equal(fvec, seq(0,4) * 1.1))

ivec <- npyLoad("ivec.npy", "integer")
print(ivec)
stopifnot(all.equal(ivec, seq(0,4)))


