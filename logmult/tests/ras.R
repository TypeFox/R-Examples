tab <- matrix(round(runif(12, 0, 2)), 3, 4)

library(logmult)
data(gss8590)

tab <- gss8590[,,1]
row <- c(1, 1, 1, 1)
col <- c(1, 1, 1, 1, 1) * 4/5
tab2 <- ras(tab, row, col)
lambda1 <- logmult:::lambda(tab)
lambda2 <- logmult:::lambda(tab2)
stopifnot(all.equal(lambda1, lambda2))
