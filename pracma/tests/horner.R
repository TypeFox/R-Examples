##
##  h or n e r . R  Test Suite
##


horner <- pracma::horner
hornerdefl <- pracma::hornerdefl

p <- c(1, 0, 1)
x <- c(-2, -1, 0, 1, 2)
identical(horner(p, x)$y, x^2 + 1)
identical(horner(p, x)$dy, 2*x)

p <- c(1, -6, 11, -6)
identical(hornerdefl(p, 3)$y, 0)
identical(hornerdefl(p, 3)$q, (c(1, -3, 2)))
