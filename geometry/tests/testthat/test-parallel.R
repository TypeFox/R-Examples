context("Interaction with parallel package")
library(parallel)
test_that("delaunayn can be called with mc.apply", {
  ## mc.cores must be 1 on Windows. Otherwise use only 2 cores to comply
  ## with CRAN guidelines.
  mc.cores <- ifelse(Sys.info()[1] == "Windows", 1, 2)

  ## Set seed for replicability
  set.seed(1)

  ## Create points and try standard Delaunay Triangulation
  N <- 100000
  P <- matrix(runif(2*N), N, 2)
  T <- delaunayn(P)
  expect_that(nrow(T), equals(199966))

  ## Now try out the parallel version. 
  Ts <- mclapply(list(P, P, P, P), delaunayn, mc.cores=mc.cores)
  expect_that(length(Ts), equals(4))
  expect_that(nrow(Ts[[1]]), equals(199966))
  expect_that(Ts[[1]], equals(T))
})
