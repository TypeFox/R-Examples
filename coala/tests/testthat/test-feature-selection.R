context("Feature Selection")

test_that("generation of initial selection cmd works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  model  <- model_theta_tau() +
    feat_selection(111, 222, 333, population = "all", time = 5)
  cmd <- msms$get_cmd(model)
  expect_true(grepl(" -N 10000 ", cmd))
  expect_true(grepl(" -SForceKeep ", cmd))
  expect_true(grepl(" -SI 5 2 5e-04 5e-04 ", cmd))
  expect_true(grepl(" -Sp 0.5 ", cmd))
  expect_true(grepl(" -SAA 111 ", cmd))
  expect_true(grepl(" -SAa 222 ", cmd))
  expect_true(grepl(" -Saa 333 ", cmd))
  expect_true(grepl(" $", cmd))

  model <- model_theta_tau() +
    feat_selection(strength_A = 123, population = 1, time = 5, Ne = 1122,
                   start_frequency = c(0.1), position = 0.3)
  cmd <- msms$get_cmd(model)
  expect_true(grepl(" -Sc 0 1 123 ", cmd))
  expect_true(grepl(" -SI 5 2 0.1 0 ", cmd))
  expect_true(grepl("-N 1122", cmd))
  expect_true(grepl(" -Sp 0.3 ", cmd))
  expect_true(grepl(" $", cmd))

  model <- model_theta_tau() +
    feat_selection(strength_A = 123, population = 1, time = 5,
                   start_frequency = c(0.1, 0.2))
  expect_true(grepl(" -SI 5 2 0.1 0.2 ", msms$get_cmd(model)))
})


test_that("generation of selection cmd with multiple demes works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  model  <- model_theta_tau() +
    feat_selection(111, 222, 333, population = 1, time = 5) +
    feat_selection(444, 555, 666, population = 2, time = 1.0, start = FALSE)
  cmd <- msms$get_cmd(model)
  expect_true(grepl(" -N 10000 ", cmd))
  expect_true(grepl(" -SForceKeep ", cmd))
  expect_true(grepl(" -SI 5 2 5e-04 0 ", cmd))
  expect_true(grepl(" -Sp 0.5 ", cmd))
  expect_true(grepl(" -Sc 0 1 111 222 333 ", cmd))
  expect_true(grepl(" -Sc 1 2 444 555 666 ", cmd))
  expect_true(grepl(" $", cmd))
})


test_that("msms can simulate selection", {
  if (!has_msms()) skip("msms not installed")
  # With one population
  model <- coal_model(5, 1, 100) + par_named("s") + par_named("t") +
    feat_selection(par_expr(2 * s), par_expr(.5 * s), time = par_expr(t * 2)) +
    feat_mutation(5) +
    sumstat_sfs()

  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model, pars = c(s = 20, t = 0.01))
  expect_that(stat$sfs, is_a("numeric"))

  # With three population
  if (!has_msms()) skip("msms not installed")
  model <- coal_model(c(5, 5, 5), 1, 100) +
    feat_selection(1000, 500, 1, time = 0.01) +
    feat_mutation(5) +
    feat_migration(1, symmetric = TRUE) +
    sumstat_sfs()
  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))

  # With additive selection
  model <- coal_model(5, 1, 100) +
    feat_selection(strength_A = 123, time = 0.03) +
    feat_mutation(5) +
    sumstat_sfs()
  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})
