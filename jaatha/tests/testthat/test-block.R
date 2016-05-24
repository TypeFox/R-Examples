context("Block")

border <- matrix(c(0.4, 0.2, 0.6, 0.4), 2, 2)
border_1dim <- matrix(c(0, 1), 1, 2)

test_that("blocks can be initalized", { #nolint
  block <- create_block(border)
  expect_equivalent(block$get_border(), border)
  
  expect_error(create_block(t(border)))
  expect_error(create_block(border[, 2:1]))
  
  expect_error(create_block(matrix(c(-1, .5, .5, 2), 2)))
  block <- create_block(matrix(c(-1, .5, .5, 2), 2), TRUE)
  expect_equivalent(block$get_border(), matrix(c(0, .5, .5, 1), 2))
})


test_that("test if a block includes a point works", {
  block <- create_block(border)
  expect_equal(block$includes(c(0.5, 0.3)), TRUE)
  expect_equal(block$includes(c(0.4, 0.3)), TRUE)
  expect_equal(block$includes(c(0.6, 0.3)), TRUE)
  expect_equal(block$includes(c(0.6, 0.4)), TRUE)
  expect_equal(block$includes(c(0.3, 0.3)), FALSE)
  expect_equal(block$includes(c(0.3, 0.1)), FALSE)
  expect_equal(block$includes(c(0.4, 0.61)), FALSE)
})


test_that("block returns its middle", {
  block <- create_block(border)
  expect_equivalent(block$get_middle(), c(.5, .3))
  
  block <- create_block(border_1dim)
  expect_equivalent(block$get_middle(), .5)
})


test_that("block returns its corners", {
  block <- create_block(border)
  expect_equivalent(block$get_corners(), matrix(c(.4, .6, .4, .6,
                                                  .2, .2, .4, .4), 4, 2))
  
  block <- create_block(border_1dim)
  expect_equivalent(block$get_corners(), matrix(c(0, 1), 1, 2))
})


test_that("sampling of parameters works", {
  block <- create_block(border)
  set.seed(17)
  random_pars <- block$sample_pars(10)
  expect_that(random_pars, is_a("matrix"))
  expect_equal(dim(random_pars), c(10, 2))
  expect_true(all(apply(random_pars, 1, function(x) block$includes(x))))
  
  set.seed(17)
  random_pars_2 <- block$sample_pars(10, TRUE)
  expect_equal(random_pars_2, rbind(random_pars, block$get_corners()))
  expect_true(all(apply(random_pars, 1, function(x) block$includes(x))))
  
  block <- create_block(border_1dim)
  random_pars <- block$sample_pars(5)
  expect_that(random_pars, is_a("matrix"))
  expect_equal(dim(random_pars), c(5, 1))
  expect_true(all(apply(random_pars, 1, function(x) block$includes(x))))
})


test_that("it removes the outer area of a block", {
  border <- matrix(c(0.4, 0.1, 0.6, 0.4), 2, 2)
  block <- create_block(border)
  
  expect_equivalent(block$get_interior(0.1),
                    matrix(c(0.42, 0.13, 0.58, 0.37), 2))
  expect_equivalent(block$get_interior(0.25),
                    matrix(c(0.45, 0.175, 0.55, 0.325), 2))
  
  border <- matrix(c(0.1, 0.9), 1, 2)
  block <- create_block(border)
  expect_equivalent(block$get_interior(0.1),
                    matrix(c(0.18, 0.82), 1))
  expect_equivalent(block$get_interior(0.2),
                    matrix(c(0.26, 0.74), 1))
})
