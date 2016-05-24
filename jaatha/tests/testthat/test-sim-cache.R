context("Simulation Cache")

sim_data <- lapply(1:10 / 10, function(x) list(pars_normal = rep(x, 3)))
block_all <- create_block(matrix(c(0, 1), 3, 2, byrow = TRUE))
block_small <- create_block(matrix(c(0, .25), 3, 2, byrow = TRUE))

test_that("caching of simulations works", {
  sc <- create_sim_cache()
  expect_equal(sc$get_size(), 0)
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 10)
  expect_equal(sc$get_sim_data(block_all), sim_data)
  expect_equal(sc$get_sim_data(block_small), sim_data[1:2])
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 20)
  expect_equal(sc$get_sim_data(block_small), sim_data[c(1:2, 1:2)])
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 30)
  expect_equal(sc$get_sim_data(block_small), sim_data[c(1:2, 1:2, 1:2)])
  
  expect_error(sc$add(5))
  capture.output(expect_error(sc$add(list(5))))
  capture.output(expect_error(sc$add(list(list(5)))))
  
  # One parameter
  sc <- create_sim_cache()
  sim_data <- lapply(1:10 / 10, function(x) list(pars_normal = x))
  sc$add(sim_data)
  expect_equal(sc$get_sim_data(create_block(matrix(0:1, 1, 2))), sim_data)
})

test_that("limiting of cache works", {
  sc <- create_sim_cache(5)
  block_all <- create_block(matrix(c(0, 1), 3, 2, byrow = TRUE))
  expect_equal(sc$get_size(), 0)
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 5)
  expect_equal(sc$get_sim_data(block_all), sim_data[1:5])
  
  sc <- create_sim_cache(15)
  sc$add(sim_data)
  expect_equal(sc$get_size(), 10)
  expect_equal(sc$get_sim_data(block_all), sim_data)
  sc$add(sim_data)
  expect_equal(sc$get_size(), 15)
  expect_equal(sc$get_sim_data(block_all), sim_data[c(1:10, 1:5)])
})

test_that("the cache can be disabled", {
  sc <- create_sim_cache(0)
  sc$add(sim_data)
  expect_equal(sc$get_sim_data(block_all), sim_data)
})
