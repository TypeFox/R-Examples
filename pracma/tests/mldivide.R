##
##  m l d i v i d e . R  test suite
##


mldivide <- pracma::mldivide
mrdivide <- pracma::mrdivide

A <- matrix(c(8,1,6, 3,5,7, 4,9,2), nrow = 3, ncol = 3, byrow = TRUE)
identical(all.equal(mldivide(A, A), diag(1, 3, 3), tolerance=1e-7), TRUE)
identical(all.equal(mrdivide(A, A), diag(1, 3, 3), tolerance=1e-7), TRUE)
