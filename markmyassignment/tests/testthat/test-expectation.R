
context("expectations")

test_that(desc="expect_function_self_contained()",{
  h <- "1"
  f <- function() h
  g <- function() h <- 1; h
  expect_failure(expect_function_self_contained(object = f))
  expect_success(expect_function_self_contained(object = g))
})


test_that(desc="expect_attached_package()",{
  expect_success(expect_attached_package("base"))
  expect_failure(expect_attached_package("advfsda"))
})


test_that(desc="expect_function_arguments()",{
  f <- function(x, y) x^2
  g <- function() 3
  expect_failure(expect_function_arguments(f, "x"))
  expect_success(expect_function_arguments(f, c("x", "y")))
  expect_failure(expect_function_arguments(f, c("x", "y", "z")))
  expect_success(expect_function_arguments(g, NULL))
})


test_that(desc="expect_function_code()",{
  expect_success(expect_function_code(object = base::mean, expected = "UseMethod"))
  expect_failure(expect_function_code(object = base::mean, expected = "markmyassignment"))
})

