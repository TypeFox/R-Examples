context("Parameter Zero Inflation")

test_that("parameters with zero inflation can be initialized", {
  set.seed(114422)
  par <- par_zero_inflation(par_const(5), 0.20)
  expect_true(is.par(par))
  expect_true(is.par_variation(par))
  expect_is(par$get_expression(), "expression")
  x <- vapply(1:10000, function(x) par$eval(), numeric(1))
  expect_that(mean(x), is_more_than(3.5))
  expect_that(mean(x), is_less_than(4.5))

  par <- par_zero_inflation(5, .2)
  expect_that(par$eval(), is_a("numeric"))
  par <- par_zero_inflation("5", .2)
  expect_that(par$eval(), is_a("numeric"))
  par <- par_zero_inflation(par_named("x"), .2)
  x <- 100
  expect_that(par$eval(), is_a("numeric"))

  par <- par_zero_inflation(5, "a")
  a <- .1
  expect_that(par$eval(), is_a("numeric"))
  par <- par_zero_inflation(5, "a")
  a <- .1
  expect_that(par$eval(), is_a("numeric"))
})


test_that("non-random zero inflation works", {
  model <- coal_model(5, 100) +
    feat_mutation(par_zero_inflation(5, 0.211, FALSE))
  expect_true(has_variation(model))
  tmplt <- scrm_create_cmd_template(model)
  cmd_tmplt <- fill_cmd_template(tmplt, model, NULL, locus_group = 1)
  expect_equivalent(cmd_tmplt, data.frame(locus_number = c(21, 79),
                                          command = c("-t 0 ", "-t 5 "),
                                          stringsAsFactors = FALSE))

  model <- coal_model(5, 100) +
    feat_mutation(par_zero_inflation(5, 0.4, FALSE))
  expect_true(has_variation(model))
  tmplt <- scrm_create_cmd_template(model)
  cmd_tmplt <- fill_cmd_template(tmplt, model, NULL, locus_group = 1)
  expect_equivalent(cmd_tmplt, data.frame(locus_number = c(40, 60),
                                          command = c("-t 0 ", "-t 5 "),
                                          stringsAsFactors = FALSE))
})


test_that("Adding a parameter with variation to a model works", {
  expect_false(has_variation(coal_model(5)))

  model <- coal_model(5) + par_zero_inflation(5, 100)
  expect_equal(length(get_parameter(model)), 0)
  expect_true(has_variation(model))

  model <- coal_model(5) + par_zero_inflation(par_named("x"), 100)
  expect_equal(length(get_parameter(model)), 1)
  expect_equal(get_parameter(model)[[1]]$get_name(), "x")
  expect_true(has_variation(model))
})
