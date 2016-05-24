##
##  ha d a m a r d . R  Test suite
##


hadamard <- pracma::hadamard
Toeplitz <- pracma::Toeplitz

all.equal(hadamard(2),
          matrix(c(1, 1, 1, -1), 2, 2))
all.equal(hadamard(4),
          matrix(c(1,  1,  1,  1,
                   1, -1,  1, -1,
                   1,  1, -1, -1,
                   1, -1, -1,  1), 4, 4))
# H12 <- hadamard(12)
# all.equal(t(H12) %*% H12,
#           diag(12, 12, 12))
# H20 <- hadamard(20)
# all.equal(t(H20) %*% H20,
#           diag(20, 20, 20))

all.equal(Toeplitz(c(1, 2, 4, 6, 8), c(1, 3, 5, 7, 9)),
          matrix(c(1, 3, 5, 7, 9,
                   2, 1, 3, 5, 7,
                   4, 2, 1, 3, 5,
                   6, 4, 2, 1, 3,
                   8, 6, 4, 2, 1), 5, 5, byrow = TRUE))
