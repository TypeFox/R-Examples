##
##  m t r a c e . R  Test suite
##

Diag <- pracma::Diag

all.equal(Diag(matrix(1:12,3,4),  1), c(4,8,12))
all.equal(Diag(matrix(1:12,3,4), -1), c(2,6))
identical(Diag(Diag(c(1,5,9)), 0), c(1,5,9))
