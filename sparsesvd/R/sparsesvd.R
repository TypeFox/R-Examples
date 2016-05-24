sparsesvd <- function (M, rank=0L, tol=1e-15, kappa=1e-6) {
  if (is.matrix(M)) stop("argument must be a sparse real matrix")
  if (!is(M, "dMatrix")) stop("only sparse real dMatrix format (from Matrix package) is currently supported")
  if (!is(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  .Call(svdLAS2_, dim(M), M@i, M@p, M@x, as.integer(rank), as.double(tol * c(-1, 1)), as.double(kappa))
}
