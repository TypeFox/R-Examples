
library("slam")

##
x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3)
x
sx <- as.simple_triplet_matrix(x)

y <- matrix(1:6, nrow = 3)
sy <- as.simple_triplet_matrix(y)

identical(tcrossprod(x, y), tcrossprod_simple_triplet_matrix( x, sy))
identical(tcrossprod(x),    tcrossprod_simple_triplet_matrix(sx))
identical(tcrossprod(x, y), tcrossprod_simple_triplet_matrix(sx, sy))
identical(tcrossprod(x, y), tcrossprod_simple_triplet_matrix(sx,  y))

identical(crossprod(x, y),  crossprod_simple_triplet_matrix( x, sy))
identical(crossprod(x),     crossprod_simple_triplet_matrix(sx))
identical(crossprod(x, y),  crossprod_simple_triplet_matrix(sx, sy))
identical(crossprod(x, y),  crossprod_simple_triplet_matrix(sx,  y))

identical(crossprod(x, y),  matprod_simple_triplet_matrix(t( x), sy))
identical(crossprod(x, y),  matprod_simple_triplet_matrix(t(sx), sy))
identical(crossprod(x, y),  matprod_simple_triplet_matrix(t(sx),  y))

## Note that correctness under bailout is covered elsewhere.

##
