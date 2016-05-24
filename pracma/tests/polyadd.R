##
##  p o l y m u l . R  Test suite
##


polyadd <- pracma::polyadd

identical(polyadd(c(1, 1, 1), 1), c(1, 1, 2))
identical(polyadd(c(1, 1, 1), c(0, 1)), c(1, 1, 2))
identical(polyadd(c(0.5, 1, 1), c(0.5, 1, -1)), c(1, 2, 0))
identical(polyadd(c(0.5, 1, 1), c(-0.5, -1, 1)), c(2))
identical(polyadd(c(0, 0, 1, 2, 2), c(0, 1, 2, 3, 4)), c(1, 3, 5, 6))
