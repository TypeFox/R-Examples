context("Feature Migration")

test_that("Creation of migration features works", {
  feat <- feat_migration(2, 2, 1)
  expect_equal(feat$get_rate(), "par(2)")
  expect_equal(feat$get_population(), c(from = 2, to = 1))
  expect_equal(feat$get_time(), "par(0)")

  feat <- feat_migration(2, symmetric = TRUE, time = 5)
  expect_equal(feat$get_population(), "all")

  expect_error(feat_migration(3, "A", 5))
  expect_error(feat_growth(3, 1:2, 5))
})


test_that("generating scrm cmd for growth works", {
  model <- coal_model(4:5, 1) + feat_migration(par_range("m", 1, 2), 2, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model),
               "scrm 9 1 -I 2 4 5 -em 0 1 2 m ")
  model <- coal_model(4:5, 1) + feat_migration(3, symmetric = TRUE, time = 5)
  expect_equal(get_simulator("scrm")$get_cmd(model),
               "scrm 9 1 -I 2 4 5 -eM 5 3 ")
})


test_that("migration can be simulated with scrm", {
  model <- coal_model(4:5, 1) +
    feat_mutation(1) +
    par_range("m", 1, 2) +
    feat_migration(par_expr(2 * m), 2, 1, time = par_expr(log(m)))
  expect_that(get_simulator("scrm")$simulate(model, c(m = 1.5)), is_a("list"))
})
