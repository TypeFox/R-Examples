##
##  c o m b s . R  Test Suite
##

combs <- pracma::combs
randcomb <- pracma::randcomb

identical(combs(2, 1), 2)
identical(combs(c(1, 2, 3), 2), matrix(rep(c(1, 2, 3), each = 2), 3, 2))
identical(nrow(combs(1:6, 4)), 15L)

all(c(1,2,3) %in% randcomb(c(1,2,3), 3))
