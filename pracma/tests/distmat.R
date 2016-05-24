##
##  d i s t m a t . R  tests
##

distmat <- pracma::distmat

A <- c(0.0, 0.0, 0.0)
B <- matrix(c(
        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        0, 1, 1,
        1, 0, 1,
        1, 1, 0,
        1, 1, 1),
     nrow = 8, ncol = 3, byrow = TRUE)

all.equal(drop(distmat(A, B)),
          c(0, 1, 1, 1, sqrt(2), sqrt(2), sqrt(2), sqrt(3)))
