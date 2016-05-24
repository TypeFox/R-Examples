##
##  p a s c a l . R  Test Suite
##

pascal <- pracma::pascal
nchoosek <- pracma::nchoosek

identical(pascal(3),
          matrix(c(1, 1, 1, 1, 2, 3, 1, 3, 6), 3, 3))

identical(nchoosek(6, 1), choose(6, 1))
identical(nchoosek(6, 2), choose(6, 2))
identical(nchoosek(6, 3), choose(6, 3))
identical(nchoosek(6, 4), choose(6, 4))
identical(nchoosek(6, 5), choose(6, 5))
identical(nchoosek(6, 6), choose(6, 6))
