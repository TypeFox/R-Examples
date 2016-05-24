context("Stat Cubes")


test_that("calculation of break_values works", {
  calc_func <- function(x) matrix(x, ncol = 1)
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE, breaks = .5)
  expect_equal(stat$generate_data_opts(1:20), list(break_values = list(10.5)))
  
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE, 
                             breaks = c(.25, .5, .75))
  expect_equal(stat$generate_data_opts(1:20), 
               list(break_values = list(c(5.75, 10.5, 15.25))))
  
  calc_func <- function(x) matrix(x, ncol = 2)
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE, breaks = .5)
  expect_equal(stat$generate_data_opts(1:20), 
               list(break_values = list(5.5, 15.5)))
  
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE, 
                             breaks = c(.25, .75))
  expect_equal(stat$generate_data_opts(1:20), 
               list(break_values = list(c(3.25, 7.75), c(13.25, 17.75))))
  
  expect_error(create_jaatha_stat("cube", calc_func, 
                                  poisson = FALSE, breaks = -1))
  expect_error(create_jaatha_stat("cube", calc_func, 
                                  poisson = FALSE, breaks = 2))
})


test_that("calculation of break_values supports vectors", {
  stat <- create_jaatha_stat("cube", I, poisson = FALSE, breaks = .5)
  expect_equal(stat$generate_data_opts(1:20), list(break_values = list(10.5)))
  stat <- create_jaatha_stat("cube", I, poisson = FALSE, 
                             breaks = c(.25, .5, .75))
  expect_equal(stat$generate_data_opts(1:20), 
               list(break_values = list(c(5.75, 10.5, 15.25))))
})


test_that("generation of cubes works", {
  value <- c(1:6, 1:6, 1:6)
  calc_func <- function(x) matrix(x, ncol = 3)
  opts <- list(break_values = list(1:2 + .5, 3.5, c(1, 3, 5) + .5))
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE)
  
  cube <- array(stat$calculate(value, opts), c(3, 2, 4))
  expect_equal(sum(cube), 6)
  expect_equal(cube[1, 1, 1], 1)
  expect_equal(cube[2, 1, 2], 1)
  expect_equal(cube[3, 1, 2], 1)
  expect_equal(cube[3, 2, 3], 2)
  expect_equal(cube[3, 2, 4], 1)
  
  
  # 2D
  value <- c(1:6, 1:6)
  calc_func <- function(x) matrix(x, ncol = 2)
  opts <- list(break_values = list(1:2 + .5, 3.5))
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE)
  
  cube <- array(stat$calculate(value, opts), c(3, 2))
  expect_equal(sum(cube), 6)
  expect_equal(cube[1, 1], 1)
  expect_equal(cube[2, 1], 1)
  expect_equal(cube[3, 1], 1)
  expect_equal(cube[3, 2], 3)
})


test_that("cubes are calculated with NaNs present", {
  calc_func <- function(x) matrix(x, ncol = 2)
  opts <- list(break_values = list(1:2 + .5, 3.5))
  stat <- create_jaatha_stat("cube", calc_func, poisson = FALSE)
  
  value <- c(1:5, NaN, 1:6)
  cube <- array(stat$calculate(value, opts), c(3, 2))
  expect_equal(sum(cube), 5)
  
  value <- c(1:5, NA, 1:6)
  cube <- array(stat$calculate(value, opts), c(3, 2))
  expect_equal(sum(cube), 5)
})
