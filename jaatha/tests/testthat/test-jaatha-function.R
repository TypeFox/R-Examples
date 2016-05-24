context("Jaatha Estimation")

test_that("the main function works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  results <- jaatha(model, data, repetitions = 2, sim = 10, cores = 1, 
                    max_steps = 15)
  
  expect_is(results, "list")
  expect_true(is.finite(results$loglikelihood))
  expect_true(all(results$param > 1))
  expect_true(is_single_logical(results$converged))
})


test_that("output can be suppessed by verbose argument", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  expect_message(
    jaatha(model, data, repetitions = 1, sim = 10, cores = 1, 
           max_steps = 5, verbose = FALSE),
    NA
  )
  
  expect_message(
    jaatha(model, data, repetitions = 1, sim = 10, cores = 1, 
           max_steps = 5, verbose = TRUE),
    "Step"
  )
})


test_that("it supports a one parameter model", {
  model <- create_jaatha_model(function(x) rpois(10, x),
                               par_ranges = matrix(c(0.1, 10), 1, 2),
                               sum_stats = list(stat_identity(), stat_sum()),
                               test = FALSE)
  
  data <- create_test_data(model)
  results <- jaatha(model, data, repetitions = 1, sim = 10, cores = 1, 
                    max_steps = 5 * 3)
  
  expect_is(results, "list")
  expect_true(is.finite(results$loglikelihood))
  expect_true(all(results$estimate > 1))
  expect_identical(results$args$model, model)
  expect_identical(results$args$data, data)
  expect_equal(results$args$repetitions, 1)
  expect_equal(results$args$sim, 10)
  expect_equal(results$args$cores, 1)
  expect_equal(results$args$max_steps, 15)
})
