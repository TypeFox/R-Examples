context("Parameter Priors")

test_that("creating parameters with priors works", {
  par <- par_prior("x", rnorm(1))
  expect_equal(par$get_name(), "x")
  expect_true(is.prior_par(par))

  set.seed(18)
  x <- sapply(1:10, function(i) par$sample())
  expect_equal(length(unique(x)), 10)

  set.seed(18)
  y <- sapply(1:10, function(i) par$sample())
  expect_equal(x, y)
})


test_that("sampling parameter priors works", {
  expect_equal(sample_par_priors(coal_model(5)), numeric(0))

  model <- coal_model(5) + par_prior("a", 1.5) + par_prior("b", 2.5)
  expect_equal(sample_par_priors(model), c(a = 1.5, b = 2.5))
})
