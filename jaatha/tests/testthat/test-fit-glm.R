context("GLM Fitting")


test_that("fit_glm works for PoiInd", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  
  glms <- fit_glm(model$get_sum_stats()[[1]], sim_data)
  expect_true(is.list(glms))
  expect_equal(10, length(glms))
  
  glms <- fit_glm(model$get_sum_stats()[[2]], sim_data)
  expect_true(is.list(glms))
  expect_equal(1, length(glms))
  
  expect_error(fit_glm(1:10, sim_data))
})


test_that("fit_glm works", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  
  glm_fit <- fit_glm(model, sim_data)
  expect_that(glm_fit, is_a("list"))
  expect_equal(length(glm_fit), 2)
  expect_equal(names(glm_fit), c("id", "sum"))
})
