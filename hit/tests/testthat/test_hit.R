context("HIT testing")

test_that("hit testing", {
  set.seed(123)
  n <- 100
  p <- 150
  # x with correlated columns
  corMat <- toeplitz((p:1 / p)^3)
  corMatQ <- chol(corMat)
  x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
  colnames(x) <- paste0("x", 1:p)
  # y
  y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
  # hierarchy
  dend <- as.dendrogram(hclust(dist(t(x))))
  hier <- as.hierarchy(dend)
  # HIT
  fit <- hit(x, y, hier)
  # checks
  expect_equal(class(fit), "hit")
  expect_equal(unname(unlist(lapply(fit, class))), 
               c("numeric", "numeric","hierarchy", 
                 "logical", "numeric", "numeric", "numeric"))
  expect_equal(names(fit), 
               c("pValues", "selectFreq", "hierarchy", 
                 "tested", "alpha", "lambda", "max.p.esti"))
})
