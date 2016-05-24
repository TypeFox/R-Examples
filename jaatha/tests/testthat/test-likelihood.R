context("Likelihood Estimation")

test_that("llh is approximatied for basic statistics", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  llh <- approximate_llh(model$get_sum_stats()[[1]], data, c(.5, .5), 
                         glms, 10, 1)
  expect_true(is.numeric(llh))
  expect_true(llh <= 0)
  
  llh2 <- approximate_llh(model$get_sum_stats()[[1]], data, c(.5, .5), 
                          glms, 10, 2)
  expect_true(is.numeric(llh2))
  expect_true(llh2 <= 0)
  expect_true(llh != llh2)
  
  expect_error(approximate_llh(1:3, data, c(.5, .5), glms, 10, 2))
})


test_that("llh is approximatied for complete models", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  llh <- approximate_llh(model, data, c(.5, .5), glms, 10)
  expect_true(is.numeric(llh))
  expect_true(llh <= 0)
  
  llh1 <- approximate_llh(model$get_sum_stats()[[1]], 
                          data, c(.5, .5), glms, 10, 1)
  llh2 <- approximate_llh(model$get_sum_stats()[[2]], 
                          data, c(.5, .5), glms, 10, 1)
  expect_equal(llh, llh1 + llh2)
})


test_that("llh optimization works", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  opt_llh <- optimize_llh(block, model, data, glms, 20)
  expect_equal(length(opt_llh$par), 2)
  expect_true(all(opt_llh$par > 0 & opt_llh$par < 1))
  expect_true(opt_llh$value < 0)
})


test_that("precise llh estimation works", {
  model <- create_test_model()
  data <- create_test_data(model)
  llh <- estimate_llh(model, data, c(.5, .5), sim = 20, 
                      cores = 1, normalized = TRUE)
  expect_equivalent(llh$param, c(.5, .5))
  expect_that(llh$value, is_less_than(0))
  
  model <- create_test_model()
  data <- create_test_data(model)
  expect_error(estimate_llh(model, data, c(1.5, 1.5), sim = 20, 
                            cores = 1, normalized = TRUE))
  llh <- estimate_llh(model, data, c(1.5, 1.5), sim = 20, 
                      cores = 1, normalized = FALSE)
  expect_equivalent(llh$param, model$get_par_ranges()$normalize(c(1.5, 1.5)))
  expect_that(llh$value, is_less_than(0))
})


test_that("it estimates local llh maxima", {
  model <- create_test_model()
  data <- create_test_data(model)
  sim_cache <- create_sim_cache()
  block <- create_block(matrix(0:1, 2, 2, byrow = TRUE))
  
  est <- estimate_local_ml(block, model, data, 20, 1, sim_cache)
  expect_is(est$par, "numeric")
  expect_equal(length(est$par), 2)
  expect_true(all(est$par >= 0 & est$par <= 1))
  
  expect_is(est$value, "numeric")
  expect_equal(length(est$value), 1)
  expect_true(all(est$value <= 0))
})
