context("Feature Size Change")

test_that("generating scrm cmd for size changes works", {
  get_cmd <- get_simulator("scrm")$get_cmd
  model <- coal_model(15, 1) + feat_size_change(par_range("alpha", 0, 1), 1)
  expect_equal(get_cmd(model), "scrm 15 1 -N alpha ")

  model <- coal_model(15, 1) + feat_size_change(5, "all")
  expect_equal(get_cmd(model), "scrm 15 1 -N 5 ")

  model <- coal_model(15, 1) + feat_size_change(5, "all", 1)
  expect_equal(get_cmd(model), "scrm 15 1 -eN 1 5 ")

  model <- coal_model(14:15, 1) + feat_size_change(5, "all")
  expect_equal(get_cmd(model), "scrm 29 1 -I 2 14 15 -N 5 ")

  model <- coal_model(14:15, 1) + feat_size_change(5, "all", 1)
  expect_equal(get_cmd(model), "scrm 29 1 -I 2 14 15 -eN 1 5 ")

  model <- coal_model(14:15, 1) + feat_size_change(5, 2)
  expect_equal(get_cmd(model), "scrm 29 1 -I 2 14 15 -n 2 5 ")

  model <- coal_model(14:15, 1) + feat_size_change(5, 2, 1)
  expect_equal(get_cmd(model), "scrm 29 1 -I 2 14 15 -en 1 2 5 ")
})


test_that("simulating a size change works", {
  scrm <- get_simulator("scrm")
  model_tmp <- coal_model(5, 1) +
    feat_mutation(1) +
    feat_size_change(2, population = 1) +
    sumstat_sfs()

  sum_stats <- scrm$simulate(model_tmp)
  expect_true(is.numeric(sum_stats$sfs))
})
