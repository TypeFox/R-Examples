##
##  p o l y m u l . R  Test suite
##


polymul <- pracma::polymul

identical(polymul(c(0.5), c(2, 4, 8)), c(1, 2, 4))
identical(polymul(c(2.5, 1.5, 0.5), c(2)), c(5, 3, 1))
identical(polymul(c(1, 1, 1), c(0, 1, 1, 1)), c(1, 2, 3, 2, 1))
identical(polymul(c(1, 0, 0), c(0, 0, 1)), c(1, 0, 0))
