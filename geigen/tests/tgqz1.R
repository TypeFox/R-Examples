
library("geigen")
source("testqz.R")

# Data from R-help mailinglist: Solve an ordinary or generalized eigenvalue problem in R?
# started on 19-04-2012

A <- matrix(c(1457.738, 1053.181, 1256.953,
              1053.181, 1213.728, 1302.838,
              1256.953, 1302.838, 1428.269), nrow=3, byrow=TRUE)

B <- matrix(c(4806.033, 1767.480, 2622.744,
              1767.480, 3353.603, 3259.680,
              2622.744, 3259.680, 3476.790), nrow=3, byrow=TRUE)

# Test interface to dgges (QZ method)

z <- gqz(A, B)
testqz(A,B,z)
z$sdim == 0

zn <- gqz(A, B, "N")
testqz(A,B,zn)
z$sdim == 0

identical(z,zn)
