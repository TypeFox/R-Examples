library(harvestr)
library(testthat)
library(plyr)
if(require(doParallel)) {
  cl=makeCluster(2)
  clusterEvalQ(cl, (library(harvestr)))
  registerDoParallel(cl)
}
{ context("parallel")
test_that('withseed is invariante under foreach', {
    seed=gather(1, seed=1234)[[1]]
    r <- foreach(seed = replicate(4, seed, simplify=F)) %do%
      withseed(seed, rnorm(1e5))
    s <- foreach(seed = replicate(4, seed, simplify=F)) %dopar%
      withseed(seed, rnorm(1e5))

    identical((r), (s)) 
    expect_true(all(laply(r, all.equal, r[[1]], check.attributes=F)))
    expect_true(all(laply(s, all.equal, s[[1]], check.attributes=F)))
})
test_that("farm is parallelizable.", {
    set.seed(123)
    seeds <- gather(4)
    a <- farm(seeds, runif(5))
    b <- farm(seeds, runif(5))
    c <- farm(seeds, runif(5), .parallel=T)
    
    expect_identical(a, b)
    expect_identical(a, c)
})
test_that('harvest is parallelizable with option', {
    seeds <- gather(100, seed=1234)
    e <- farm(seeds, rnorm(10))
    x <- harvest(e, sample, replace=T)
    z <- harvest(e, sample, replace=T, .parallel=T)
    expect_equivalent(noattr(x),noattr(z))
})
}
