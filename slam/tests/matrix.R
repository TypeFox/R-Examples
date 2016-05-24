
library("slam")

s <- simple_triplet_diag_matrix(1, nrow = 3)

identical(as.matrix(s) * Inf, as.matrix(s * Inf))
identical(as.matrix(s) *  NA, as.matrix(s *  NA_real_))

identical(as.matrix(s) * c(Inf, NA, 0), as.matrix(s * c(Inf, NA, 0)))

x1 <- matrix(c(1, Inf, 0, 1), nrow = 2)
x2 <- matrix(c(1, 0,  NA, 1), nrow = 2)

identical(x1 * x2, as.matrix(as.simple_triplet_matrix(x1) * x2))
identical(x1 * x2, as.matrix(as.simple_triplet_matrix(x1) *
                             as.simple_triplet_matrix(x2)))

x <- matrix(1, nrow = 3, ncol = 3)
identical(x * as.matrix(s), as.matrix(s * as.simple_triplet_matrix(x)))
identical(x / as.matrix(s), as.matrix(as.simple_triplet_matrix(x) / s))

identical(x * as.matrix(s), as.matrix(s * x))
identical(x / as.matrix(s), as.matrix(x / s))

###
