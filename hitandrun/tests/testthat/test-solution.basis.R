context("solution.basis")

# A 3-dimensional original space
n <- 3

test_that("reducing dimension by 1 works", {
  # x_1 + x_2 + x_3 = 1
  eq.constr <- 
    list(constr=t(rep(1, n)), dir='=', rhs=1)

  basis <- solution.basis(eq.constr)

  expect_equal(sum(basis$translate), 1)
  expect_equal(nrow(basis$basis), 3)
  expect_equal(ncol(basis$basis), 2) # Dimension reduced to 2

  y <- rbind(rnorm(100, 0, 100), rnorm(100, 0, 100))
  x <- basis$basis %*% y + basis$translate
  expect_equal(apply(x, 2, sum), rep(1, 100))
})

test_that("reducing dimension by 2 works", {
  # 2 x_2 = x_1
  eq.constr <- mergeConstraints(
    list(constr=rep(1, n), dir='=', rhs=1),
    list(constr=c(-1, 2, 0), dir='=', rhs=0))

  basis <- solution.basis(eq.constr)
  expect_equal(sum(basis$translate), 1)
  expect_equal(basis$translate[1], 2 * basis$translate[2])
  expect_equal(nrow(basis$basis), 3)
  expect_equal(ncol(basis$basis), 1) # Dimension reduced to 1

  y <- t(rnorm(100, 0, 100))
  x <- basis$basis %*% y + basis$translate
  expect_equal(apply(x, 2, sum), rep(1, 100))
  expect_equal(x[1,], 2 * x[2,])
})

test_that("reducing dimension by 3 works", {
  # 2 x_3 = x_2
  eq.constr <- mergeConstraints(
    list(constr=rep(1, n), dir='=', rhs=1),
    list(constr=c(-1, 2, 0), dir='=', rhs=0),
    list(constr = c(0, -1, 2), dir='=', rhs=0))

  basis <- solution.basis(eq.constr)
  expect_equal(sum(basis$translate), 1)
  expect_equal(basis$translate[1], 2 * basis$translate[2])
  expect_equal(basis$translate[2], 2 * basis$translate[3])
  expect_equal(nrow(basis$basis), 3)
  expect_equal(ncol(basis$basis), 0) # Dimension reduced to 0
})
