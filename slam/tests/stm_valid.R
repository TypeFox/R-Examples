
library("slam")
set.seed(20110217)

###
x <- matrix(sample(c(0,1), 12, TRUE), ncol = 3L)
s <- as.simple_triplet_matrix(x)
s

## make invalid row indexes
s$i[sample(seq_along(s$i), 3)] <- 0L

try(row_sums(s), silent = FALSE)

###

