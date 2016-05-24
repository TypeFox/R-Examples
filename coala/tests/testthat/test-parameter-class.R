context("Parameter Class")

test_that("Getting and Setting Expressions works", {
  expect_error(parameter_class$new(2 * x))
  expect_error(parameter_class$new(2))
  expect_error(parameter_class$new("2"))
  basic_par <- parameter_class$new(expression(2 * x))
  x <- 5
  expect_equal(basic_par$eval(), 10)

  test_env <- new.env()
  test_env[["x"]] <- 6
  expect_equal(basic_par$eval(envir = test_env), 12)
  expect_equal(basic_par$eval(), 10)

  expr <- basic_par$get_expression()
  expect_is(expr, "expression")
  expect_equal(eval(expr), 10)

  expect_true(is.par(basic_par))
  expect_false(is.named_par(basic_par))
})


test_that("par_expr works", {
  basic_par <- par_expr(2 * x)
  expect_true(is.par(basic_par))
  x <- 5
  expect_equal(basic_par$eval(), 10)
  x <- 6
  expect_equal(basic_par$eval(), 12)

  basic_par <- par_expr(sqrt(y))
  y <- 4
  expect_equal(basic_par$eval(), 2)
})


test_that("par_const works", {
  x <- 5
  basic_par <- par_const(2 * x)
  x <- 6
  expect_equal(basic_par$eval(), 10)
})


test_that("par_named works", {
  par <- par_named("theta")
  expect_true(is.par(par))
  expect_true(is.named_par(par))

  expect_equal(par$get_name(), "theta")
  theta <- 5
  expect_equal(par$eval(), 5)

  expect_equal(5, par$generate_value(c(theta = 5)))
  expect_equal(5, par$generate_value(c(x = 2, theta = 5)))
  expect_error(par$generate_value(5))

  expect_true(par$check_value(1))
  expect_true(par$check_value(2))
})


test_that("par_range works", {
  par <- par_range("theta", 1, 2)
  expect_true(is.par(par))
  expect_true(is.named_par(par))

  expect_equal(par$get_name(), "theta")
  theta <- 5
  expect_equal(par$eval(), 5)

  expect_equal(par$get_range(), 1:2)

  expect_true(par$check_value(1))
  expect_true(par$check_value(1.4))
  expect_true(par$check_value(2))
  expect_true(par$check_value(1 - 1e-10))
  expect_true(par$check_value(2 + 1e-10))

  expect_error(par$check_value("1"))
  expect_error(par$check_value(0))
  expect_error(par$check_value(3))
  expect_error(par$check_value(1:2))

  expect_error(par_range("theta", 1:2))
  expect_error(par_range("theta", 2, 1))
})


test_that("Adding an expression par to a model throws no error", {
  coal_model(5:6, 10, 100) + par_expr(2 * theta)
  coal_model(5:6, 10, 100) + par_expr(2 * theta) + par_expr(5)
})


test_that("Creation of parameter enviroment works", {
  # With named parameters
  model <- coal_model(5) + par_named("x")
  par_envir <- create_par_env(model, c(x = 5))
  expect_equal(par_envir[["x"]], 5)
  par_envir <- create_par_env(model, c(y = 2, x = 5))
  expect_equal(par_envir[["x"]], 5)
  expect_error(create_par_env(model, numeric()))
  expect_error(create_par_env(model, 1:2))
  expect_error(create_par_env(model, c(y = 2)))

  # Without parameters
  par_envir <- create_par_env(coal_model(5), numeric(0))

  # With ranged parameters (not really needed)
  par_envir <- create_par_env(model_theta_tau(), c(tau = 1, theta = 5))
  expect_equal(par_envir[["tau"]], 1)
  expect_equal(par_envir[["theta"]], 5)

  par_envir <- create_par_env(model_theta_tau(), c(theta = 5, tau = 1))
  expect_equal(par_envir[["tau"]], 1)
  expect_equal(par_envir[["theta"]], 5)

  # Additional options
  par_envir <- create_par_env(model_theta_tau(), c(tau = 1, theta = 5),
                              locus = 17)
  expect_equal(par_envir[["locus"]], 17)

  par_envir <- create_par_env(model_theta_tau(), c(tau = 1, theta = 5),
                              locus = 23, seed = 115)
  expect_equal(par_envir[["locus"]], 23)
  expect_equal(par_envir[["seed"]], 115)


  # For cmd printing
  par_envir <- create_par_env(model_theta_tau(), for_cmd = TRUE)
})


test_that("preparing parameters works", {
  expect_equal(prepare_pars(numeric(0), coal_model(5)), numeric(0))
  model <- coal_model(5) + par_named("x")
  expect_equal(prepare_pars(1.25, model), c(x = 1.25))

  model <- coal_model(5, 1) + par_prior("r", stats::rbinom(1, 3, .5))
  expect_error(prepare_pars(c("1", "2"), model))
  pars <- prepare_pars(numeric(), model)
  expect_equal(names(pars), "r")
  expect_true(all(pars %in% 0:3))

  model <- coal_model(5, 1) +
    par_prior("m", stats::rbinom(1, 3, .5)) +
    par_prior("r", stats::rbinom(1, 3, .5))
  pars <- prepare_pars(numeric(), model)
  expect_equal(names(pars), c("m", "r"))
  expect_true(all(pars %in% 0:3))

  model <- coal_model(5, 1) +
    par_named("m") +
    par_prior("r", stats::rbinom(1, 3, .5))
  pars <- prepare_pars(1, model)
  expect_equal(names(pars), c("m", "r"))
  expect_true(all(pars %in% 0:3))
  pars <- prepare_pars(c(m = 1), model)
  expect_equal(names(pars), c("m", "r"))
  expect_true(all(pars %in% 0:3))
})
