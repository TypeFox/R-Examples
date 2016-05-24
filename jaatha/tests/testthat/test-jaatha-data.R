context("Jaatha Data")


test_that("default creation of data works", {
  model <- create_test_model()
  real_data <- rep(c(17, 5), 5)
  
  jaatha_data <- create_jaatha_data(real_data, model)
  expect_true(is_jaatha_data(jaatha_data))
  
  expect_equal(jaatha_data$get_values(), 
               list(id = real_data, sum = sum(real_data)))
  expect_equal(jaatha_data$get_values("id"), real_data)
  expect_equal(jaatha_data$get_values(stat_identity()), real_data)
  
  expect_equal(jaatha_data$get_options(), list())
  expect_equal(jaatha_data$get_options("id"), NULL)
  expect_equal(jaatha_data$get_options(stat_identity()), NULL)
  
  expect_warning(jaatha_data$set_options(list(id = 5)))
  expect_equal(jaatha_data$get_options("id"), 5)
  
  log_facs_id <- log(factorial(real_data))
  log_facs_sum <- log(factorial(sum(real_data)))
  expect_equal(jaatha_data$get_log_factorial(), 
               list(id = log_facs_id, sum = log_facs_sum))
  expect_equal(jaatha_data$get_log_factorial("id"), log_facs_id)
  expect_equal(jaatha_data$get_log_factorial(stat_identity()), log_facs_id)
})


test_that("it calculates logfactorials even for large numbers", {
  model <- create_test_model()
  real_data <- c(1e6, 1:9)
  jaatha_data <- create_jaatha_data(real_data, model)
  expect_true(all(is.finite(jaatha_data$get_log_factorial("id"))))
  expect_true(all(is.finite(jaatha_data$get_log_factorial("sum"))))
})
