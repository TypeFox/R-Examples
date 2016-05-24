library(grpSLOPE)

context("grpSLOPE()")

A.vec <- c(0.26, 0.37,0.57, 0.90, 0.20, 0.89, 0.94, 0.66, 0.62,0.06)
eps   <- c(-0.7473417,-0.8858673,-0.9273622,-0.7895264,0.7119688,-0.8379027,
           -0.3327135,1.0339414,-1.2187906,-0.8921233)
A   <- diag(A.vec)
grp <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
b   <- c(0, 0, 50, 10, 0, 0, 0, 10, 0, 30)
y   <- A %*% b + eps
fdr <- 0.1
sol.c <- c(0, 0, 17.520987, 5.045372, 0, 0, 0, 0, 0, 0)
sol.beta <- c(0, 0,18.085076, 5.076808, 0, 0, 0, 0, 0, 0)
sol.group.norms <- c(0, 18.23296, 0, 0)

test_that("when the groups are consequtive blocks", {
  result <- grpSLOPE(X=A, y=y, group=grp, fdr=fdr, lambda="chiMean")
  expect_equal(result$beta, sol.beta, tolerance=1e-4)
  expect_equal(result$c, sol.c, tolerance=1e-4)
  expect_equal(as.numeric(result$group.norms), sol.group.norms, tolerance=1e-4)
  expect_identical(as.numeric(result$selected), c(1))
  expect_true(result$optimal)
  expect_is(result$iter, "numeric")
  expect_is(result$lambda, "numeric")
  expect_true(result$lambda.method == "chiMean")
  expect_is(result$sigma, "numeric")
})

test_that("when the groups are not consequtive blocks", {
  ord <- sample(1:10, 10)
  result <- grpSLOPE(X=A[ , ord], y=y, group=grp[ord], fdr=fdr, lambda="chiMean")
  expect_equal(result$beta, sol.beta[ord], tolerance=1e-4)
  expect_equal(as.numeric(result$group.norms["1"]), sol.group.norms[2], tolerance=1e-4)
  expect_true(result$optimal)
  expect_is(result$iter, "numeric")
  expect_is(result$lambda, "numeric")
  expect_true(result$lambda.method == "chiMean")
  expect_is(result$sigma, "numeric")
})

test_that("when each group is a singleton", {
  # compare to results obtained with package SLOPE
  grp <- 1:10
  sol.beta <-c(0, 0, 23.171558, 4.802997, 0, 0, 0, 4.256050, 0, 0)
  result <- grpSLOPE(X=A, y=y, group=grp, fdr=fdr, lambda="gaussian", sigma=1)
  expect_equal(result$beta, sol.beta, tolerance=1e-4)
  expect_equal(as.numeric(result$group.norms), sol.beta, tolerance=1e-4)
  expect_identical(as.numeric(result$selected), c(3, 4, 8))
  expect_true(result$optimal)
  expect_is(result$iter, "numeric")
  expect_equal(result$lambda, rep(2.575829,10), tolerance=1e-4)
  expect_true(result$lambda.method == "gaussian")
  expect_identical(result$sigma, 1)
})
