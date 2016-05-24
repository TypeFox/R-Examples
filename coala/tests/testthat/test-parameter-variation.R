context("Parameter Variation")

test_that("Parameters with variation can be initialized", {
  set.seed(114422)
  par <- par_variation(par_const(5), 100)
  expect_true(is.par(par))
  expect_true(is.par_variation(par))
  expect_is(par$get_expression(), "expression")
  x <- vapply(1:10000, function(x) par$eval(), numeric(1))
  expect_that(mean(x), is_more_than(4.5))
  expect_that(mean(x), is_less_than(5.5))
  expect_that(var(x), is_more_than(90))
  expect_that(var(x), is_less_than(110))

  par <- par_variation(5, 100)
  expect_that(par$eval(), is_a("numeric"))
  par <- par_variation("5", 100)
  expect_that(par$eval(), is_a("numeric"))
  par <- par_variation(par_named("x"), 100)
  x <- 100
  expect_that(par$eval(), is_a("numeric"))

  par <- par_variation(5, par_named("v"))
  v <- 20
  expect_that(par$eval(), is_a("numeric"))

  expect_error(par_variation(mean, 100))
})


test_that("Adding a parameter with variation to a model works", {
  expect_false(has_variation(coal_model(5)))

  model <- coal_model(5) + par_variation(5, 100)
  expect_equal(length(get_parameter(model)), 0)
  expect_true(has_variation(model))

  model <- coal_model(5) + par_variation(par_named("x"), 100)
  expect_equal(length(get_parameter(model)), 1)
  expect_equal(get_parameter(model)[[1]]$get_name(), "x")
  expect_true(has_variation(model))

  model <- coal_model(5) + par_variation(par_named("x"), par_named("y"))
  expect_equal(length(get_parameter(model)), 2)
  expect_equal(get_parameter(model)[[1]]$get_name(), "x")
  expect_equal(get_parameter(model)[[2]]$get_name(), "y")
  expect_true(has_variation(model))
})
