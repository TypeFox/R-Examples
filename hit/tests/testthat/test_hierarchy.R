context("hirarchy testing")

test_that("hirarchy testing", {
  set.seed(123)
  n <- 100
  p <- 150
  # x with correlated columns
  corMat <- toeplitz((p:1/p)^3)
  corMatQ <- chol(corMat)
  x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
  colnames(x) <- paste0("x", 1:p)
  # y
  y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
  # hierarchy
  hc <- hclust(dist(t(x)))
  hier <- as.hierarchy(hc)
  # check:
  expect_equal(class(hier), "hierarchy")
  expect_equal(typeof(hier), "list")
})

