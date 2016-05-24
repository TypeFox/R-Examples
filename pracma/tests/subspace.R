##
##  s u b s p a c e . r  Test suite
##


orth <- pracma::orth
subspace <- pracma::subspace

is.null(orth(c()))
M <- matrix(1:12, 3, 4)
all.equal(orth(M),
          matrix(c(-0.504533, -0.760776,
                   -0.574516, -0.057141,
                   -0.644497,  0.646495), 3, 2, byrow = TRUE),
          tolerance = 1e-5)

H <- pracma::hadamard(8)
A <- H[, 2:4]
B <- H[, 5:8]
all.equal(subspace(A, B), pi/2, tolerance = 1e-10)  
