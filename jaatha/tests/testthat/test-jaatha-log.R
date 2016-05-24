context("Logging")

test_that("log initialization works", {
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 100)
  expect_equivalent(log$get_estimates(1), 
                    data.frame(rep(NA, 100), NA, NA, NA, NA))
  expect_equivalent(log$get_estimates(2), 
                    data.frame(rep(NA, 100), NA, NA, NA, NA))
})


test_that("logging estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 100)
  
  estimate <- list(par = 1:model$get_par_number(), value = 0.1)
  log$log_estimate(1, 1, estimate)

  expect_equivalent(log$get_estimates(1)[1, ], 
                    data.frame(1, 1, estimate$value, 1, 2))
  
  estimate$value <- 0.2
  log$log_estimate(1, 2, estimate)
  expect_equivalent(log$get_estimates(1)[2, ], 
                    data.frame(1, 2, estimate$value, 1, 2))
  
  log$log_estimate(2, 1, estimate)
  expect_equivalent(log$get_estimates(2)[1, ], 
                    data.frame(2, 1, estimate$value, 1, 2))
})


test_that("getting best estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 100)
  for (i in c(4, 3, 1, 5, 10, 2, 7, 9, 8) / 10) {
    log$log_estimate(1, i * 10, list(par = rep(i, model$get_par_number()), 
                                     value = i))
  }
  
  for (i in 11:20 / 10) {
    log$log_estimate(2, i * 10, list(par = rep(i, model$get_par_number()), 
                                     value = -i))
  }
  
  expect_equivalent(log$get_best_estimates(1),
                    data.frame(rep = c(1, 2),
                               step = c(10, 11),
                               llh = c(1.0, -1.1),
                               p1 = c(1.0, 1.1),
                               p2 = c(1.0, 1.1)))
  
  expect_equivalent(log$get_best_estimates(2),
                    data.frame(rep = c(1, 1, 2, 2),
                               step = c(10, 9, 11, 12),
                               llh = c(1.0, 0.9, -1.1, -1.2),
                               p1 = c(1.0, 0.9, 1.1, 1.2),
                               p2 = c(1.0, 0.9, 1.1, 1.2)))
  
  # Final
  for (i in c(4, 3, 1, 5, 10, 2, 7, 9, 8) / 10) {
    log$log_estimate("final", i * 10, list(par = rep(i, model$get_par_number()),
                                           value = -i))
  }
  
  expect_equivalent(log$get_best_estimates(2, final = TRUE),
                    data.frame(rep = c(0, 0),
                               step = c(1, 2),
                               llh = c(-0.1, -0.2),
                               p1 = c(0.1, 0.2),
                               p2 = c(0.1, 0.2)))
})


test_that("it creates the results correctly", {
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 123)
  log$log_estimate("final", 1, list(par = rep(.1, model$get_par_number()),
                                    value = -1))
  log$log_estimate("final", 2, list(par = rep(.2, model$get_par_number()),
                                    value = -2))
  
  results <- log$create_results()
  expect_true(is_jaatha_result(results))
  expect_output(print(results), "p1")
  expect_equivalent(results$estimate, 
                    model$get_par_ranges()$denormalize(c(.1, .1)))
  expect_equal(results$loglikelihood, -1)
  expect_equal(results$converged, FALSE)
  expect_is(results$args, "list")
  
  log$log_convergence(1)
  expect_equal(log$create_results()$converged, FALSE)
  log$log_convergence(2)
  expect_equal(log$create_results()$converged, TRUE)
})


test_that("the output can be suppressed", {
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 123, verbose = FALSE)
  expect_message(log$log_convergence(1), NA)
  expect_message(log$log_initialization("test"), NA)
  expect_message(log$log_new_rep(1, c(.5, .5)), NA)
  expect_message(
    log$log_estimate("final", 1, list(par = rep(.1, model$get_par_number()),
                                      value = -1)),
    NA
  )
  
  model <- create_test_model()
  log <- create_jaatha_log(model, NULL, 2, 123, verbose = TRUE)
  expect_message(log$log_convergence(1))
  expect_message(log$log_initialization("test"))
  expect_message(log$log_new_rep(1, c(.5, .5)))
  expect_message(
    log$log_estimate("final", 1, list(par = rep(.1, model$get_par_number()),
                                      value = -1))
  )
})
