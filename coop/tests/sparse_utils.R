library(coop)

set.seed(1234)

m <- 30
n <- 10
prop <- .05
x <- coop:::dense_stored_sparse_mat(m, n, prop)

all.equalish <- function(x, y, fudge=.1) abs(x-y)<=fudge

stopifnot(all.equalish(sparsity(x), 1-prop))
