context("Basic GLPK")

test_that("GLPK example 1", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max)
  expect_equal( s1$optimum, 70.5, tolerance=1e-5 )
  expect_equal( s1$solution, matrix(c(0, 10.2, 9.9)), tolerance=1e-5 )
  expect_equal( s1$status, 0 )
})

test_that("GLPK example 2", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max)
  expect_equal( s1$optimum, 62.62, tolerance=1e-5 )
  expect_equal( s1$solution, matrix(c(0, 14.14, 2.02)), tolerance=1e-5 )
  expect_equal( s1$status, 0 )
})

context("Multi optimization RHS")

test_that("changing rhs for one multi problem gives correct answer", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE
  mrhs_i   <- c(1, 3)
  mrhs_val <- matrix(c(59.4, 59.4))

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mrhs_i = mrhs_i, mrhs_val = mrhs_val)
  expect_equal( s1$optimum, 61.38, tolerance=1e-5 )
  expect_equal( s1$solution, matrix(c(0, 13.86, 1.98)), tolerance=1e-5 )
  expect_equal( s1$status, 0 )
})

test_that("two rhs for multi problem gives correct answer", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(1, 30, 1)
  max <- TRUE
  mrhs_i   <- c(1, 3)
  mrhs_val <- matrix(c(60.6, 60.6, 60, 60), nrow = 2)

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mrhs_i = mrhs_i, mrhs_val = mrhs_val)
  expect_equal( s1$optimum, c(62.62, 62), tolerance=1e-5 )
  expect_equal( s1$solution, cbind(c(0, 14.14, 2.02), c(0, 14, 2)), tolerance=1e-5 )
  expect_equal( s1$status, c(0, 0) )
})

context("Multi optimization objectives")

test_that("changing obj for one multi problem gives correct answer", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE
  mobj_i   <- c(1, 2)
  mobj_val <- matrix(c(4, 2.5))
  
  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mobj_i = mobj_i, mobj_val = mobj_val)
  expect_equal( s1$optimum, 63.12, tolerance=1e-5 )
  expect_equal( s1$solution, matrix(c(11.88, 6.24, 0.0)), tolerance=1e-5 )
  expect_equal( s1$status, 0 )
})

test_that("multiple objectives work", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE
  mobj_i   <- c(2, 3)
  mobj_val <- matrix(c(7.5, 3, 4, 5), nrow=2)
  
  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mobj_i = mobj_i, mobj_val = mobj_val)
  expect_equal( s1$optimum, c(113.625, 90.3), tolerance=1e-5 )
  expect_equal( s1$solution, cbind(c(0, 15.15, 0), c(0, 10.2, 9.9)), tolerance=1e-5 )
  expect_equal( s1$status, c(0, 0) )
})

context("Multi optimization constraints")

test_that("multi problem gives correct answer for 1 problem", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(0, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE
  mmat_i <- matrix(c(3, 3, 1, 1), ncol = 2, byrow=T)
  mmat_val <- matrix(c(2, 3))

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mmat_i = mmat_i, mmat_val = mmat_val)
  expect_equal( s1$optimum, 70.5, tolerance=1e-5 )
  expect_equal( s1$solution, matrix(c(0, 10.2, 9.9)), tolerance=1e-5 )
  expect_equal( s1$status, 0 )
})

test_that("multi problem gives correct answer for 3 problems", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(0, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(60.6, 30, 60.6)
  max <- TRUE
  mmat_i <- matrix(c(3, 3, 1, 1), ncol = 2, byrow=T)
  mmat_val <- cbind(c(2, 3), c(7.5, 3.6), c(12, -1))

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mmat_i = mmat_i, mmat_val = mmat_val)
  expect_equal( s1$optimum, c(70.5, 63.125, 80.4 ), tolerance=1e-5 )
  expect_equal( s1$solution, cbind(c(0, 10.2, 9.9), c(0, 13.8875, 2.5250), c(6.6, 16.8, 0)), tolerance=1e-5 )
  expect_equal( s1$status, c(0, 0, 0) )
})

context("Full multi optimization")

test_that("rhs and constraint changing works for 3 problems", {
  obj <- c(2, 4, 3)
  mat <- matrix(c(2, 2, 1, 4, 1, 3, 2, 2, 9), nrow = 3)
  dir <- c("<=", "<=", "<=")
  rhs <- c(0, 30, 0)
  max <- TRUE
  mmat_i   <- matrix(c(3, 2, 1, 1), ncol = 2, byrow=T)
  mmat_val <- cbind(c(2, 3), c(9, 6), c(4.5, 3))
  mrhs_i   <- c(1, 3)
  mrhs_val <- matrix(c(60.6, 60.6, 55.5, 58.5, 55.5, 57), nrow = 2)

  s1 = multi_glpk_solve_LP(obj, mat, dir, rhs, max = max, mrhs_i = mrhs_i, mrhs_val = mrhs_val, mmat_i = mmat_i, mmat_val = mmat_val)
  expect_equal( s1$optimum, c(64.3875, 34.26, 53.21053), tolerance=1e-5 )
  expect_equal( s1$solution, cbind(c(0.0, 13.25625, 3.78750), c(5.31, 5.91, 0.00), c(2.289474, 12.157895, 0.0)), tolerance=1e-5 )
  expect_equal( s1$status, c(0, 0, 0) )
})

