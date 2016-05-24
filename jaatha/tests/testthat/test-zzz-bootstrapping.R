context("Bootstrapping")

test_that("the front-end for boot works", {
  skip_if_not_installed("boot")
  if (R.Version()$major == 3 && R.Version()$minor < 2.2) {
    skip("Insufficient R Version")
  }
  
  model <- create_test_model()
  data <- create_test_data(model)
  
  results <- list(estimate = c(5.05, 5.05),
                  loglikelihood = -5,
                  converged = TRUE,
                  args = list(model = model,
                              data = data,
                              repetitions = 1,
                              sim = 10,
                              max_steps = 10,
                              init_method = "middle",
                              sim_cache_limit = 0,
                              cores = 1))
  
  boot_values <- boot_jaatha(results, 10)
  expect_is(boot_values, "boot")
  expect_true(all(is.finite(boot_values$t)))
})
