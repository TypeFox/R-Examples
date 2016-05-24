context("tests on inputs")

test_that("tests for grp variable",{
  set.seed(1234); par(mar=c(0,0,0,0))
  x <- rnorm(10, mean=rep(1:5, each=2), sd=0.4)
  y <- rnorm(10, mean=rep(c(1,2), each=5), sd=0.4)
  dataFrame <- data.frame(x=x, y=y, row.names=c(1:10))
  distxy <- dist(dataFrame)
  hc <- hclust(distxy)
  
  expect_that(sort_average_r(hc),throws_error("d variable must be a dendrogram"))
})