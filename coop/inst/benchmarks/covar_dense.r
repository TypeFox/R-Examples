library(coop)
library(rbenchmark)
cols <- cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 10000
n <- 250
x <- matrix(rnorm(m*n), m, n)

benchmark(cov(x), covar(x), replications=reps, columns=cols)
