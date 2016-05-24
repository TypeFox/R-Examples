context("Extreme distances")

test_that("Main example - numeric", {
  # default distance
  mu <- c(1/7,2/7,4/7)
  nu <- c(1/4,1/4,1/2)
  x <- edistances(mu, nu)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c(0.142857142857143, 0, 0.107142857142857, 0, 0, 0.25,
                                        0, 0.285714285714286, 0.214285714285714), .Dim = c(3L, 3L), .Dimnames = list(
                                          c("1", "2", "3"), c("1", "2", "3"))))
  expect_equal(x$distances[[1]], 0.642857142857143)
  # matrix distance - named 1:3
  M <- matrix(1, nrow=3, ncol=3) - diag(3)
  rownames(M) <- colnames(M) <- 1:3
  x <- edistances(mu, nu, dist=M)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c(0.142857142857143, 0, 0.107142857142857, 0, 0, 0.25,
                                            0, 0.285714285714286, 0.214285714285714), .Dim = c(3L, 3L), .Dimnames = list(
                                              c("1", "2", "3"), c("1", "2", "3"))))
  expect_equal(x$distances[[1]], 0.642857142857143)
  # matrix distance - unnamed
  M <- matrix(1, nrow=3, ncol=3) - diag(3)
  x <- edistances(mu, nu, dist=M)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c(0.142857142857143, 0, 0.107142857142857, 0, 0, 0.25,
                                            0, 0.285714285714286, 0.214285714285714), .Dim = c(3L, 3L), .Dimnames = list(
                                              c("1", "2", "3"), c("1", "2", "3"))))
  expect_equal(x$distances[[1]], 0.642857142857143)
  # matrix distance - wrong names
  M <- matrix(1, nrow=3, ncol=3) - diag(3)
  rownames(M) <- colnames(M) <- c("a", "b", "c")
  expect_error(edistances(mu, nu, dist=M))
  # matrix distance - character
  M <- matrix("1", nrow=3, ncol=3); diag(M) <- "0"
  expect_error(edistances(mu, nu, dist=M))
})


test_that("Main example - bigq", {
  library(gmp)
  mu <- as.bigq(c(1,2,4),7)
  nu <- as.bigq(c(1,1,1),c(4,4,2))
  # default distance
  x <- edistances(mu, nu)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c("1/7", "0", "3/28", "0", "0", "1/4", "0", "2/7",
                                            "3/14"), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"),
                                                                                        c("1", "2", "3"))))
  expect_equal(x$distances[[1]], as.bigq(9,14))
  # matrix distance - numeric - unnamed
  M <- matrix(1, nrow=3, ncol=3) - diag(3)
  x <- edistances(mu, nu, dist=M)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c("1/7", "0", "3/28", "0", "0", "1/4", "0", "2/7",
                                            "3/14"), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"),
                                                                                        c("1", "2", "3"))))
  # matrix distance - character - unnamed
  M <- matrix("1", nrow=3, ncol=3); diag(M) <- "0"
  x <- edistances(mu, nu, dist=M)
  expect_true(length(x[[1]])==15 && length(x[[2]])==15)
  expect_equal(x$joinings[[1]], structure(c("1/7", "0", "3/28", "0", "0", "1/4", "0", "2/7",
                                            "3/14"), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"),
                                                                                        c("1", "2", "3"))))
  # matrix distance - bigq
  M <- as.bigq(M)
  expect_error(edistances(mu, nu, dist=M))
})
