context("Simulator scrm")

test_that("scrm can simulate seg. sites", {
  scrm <- get_simulator("scrm")

  # Generating Seg. Sites
  model <- coal_model(10, 2, 100) + feat_mutation(5) + sumstat_seg_sites()
  set.seed(110); stats_1 <- scrm$simulate(model)
  set.seed(110); stats_2 <- scrm$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$seg_sites, is_a("list"))
  expect_equal(length(stats_1$seg_sites), 2)
  expect_equal(length(get_positions(stats_1$seg_sites[[1]])),
               ncol(stats_1$seg_sites[[1]]))
})


test_that("scrm can simulate trees", {
  scrm <- get_simulator("scrm")

  model <- coal_model(10, 2, 100) + sumstat_trees()
  set.seed(110); stats_1 <- scrm$simulate(model)
  set.seed(110); stats_2 <- scrm$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
  for (i in 1:2) expect_equal(length(stats_1$trees[[i]]), 1)

  # With recombination
  set.seed(11011)
  stats <- scrm$simulate(model + feat_recombination(5))
  expect_equal(length(stats_1$trees), 2)
  for (i in 1:2) expect_true(length(stats$trees[[i]]) > 1)

  # With inter locus variation
  model <- coal_model(10, 2, 100) +
    feat_recombination(par_variation(5, 1)) +
    sumstat_trees()
  stats <- scrm$simulate(model)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
  expect_true(all(stats_1$trees[[1]] != stats_1$trees[[2]]))
})


test_that("simulation multiple loci works", {
  tmp_dir <- tempfile("scrm_test_multiple_loci")
  dir.create(tmp_dir)
  scrm <- get_simulator("scrm")
  model <- model_theta_tau() +
    feat_recombination(par_zero_inflation(1, .5)) +
    locus_averaged(10, 100) +
    sumstat_seg_sites() +
    sumstat_file(tmp_dir)
  stats <- scrm$simulate(model, c(tau = 1, theta = 5))
  expect_equal(length(stats$seg_sites), 20)
  expect_true(!any(vapply(stats$seg_sites, is.null, logical(1))))

  expect_equal(length(stats$file), 4)
  unlink(tmp_dir, recursive = TRUE)
})


test_that("printing command works", {
  scrm <- get_simulator("scrm")
  model <- model_theta_tau()
  expect_equal(scrm$get_cmd(model),
               "scrm 25 10 -I 2 10 15 -ej tau 2 1 -t theta ")
})


test_that("simulating files works", {
  scrm <- get_simulator("scrm")
  folder <- tempfile("scrm_test")
  model <- model_theta_tau() + sumstat_file(folder)
  sum_stats <- scrm$simulate(model, c(tau = 1, theta = 5))
  expect_true(is.list(sum_stats))
  expect_false(is.null(sum_stats$file))
  expect_true(is.character(sum_stats$file))
  expect_true(file.exists(sum_stats$file))
  unlink(sum_stats$file)
  unlink(folder, recursive = TRUE)
})


test_that("simulating unphased data works", {
  scrm <- get_simulator("scrm")
  model <- coal_model(5, 2, ploidy = 2) +
    feat_unphased(1) +
    feat_mutation(4) +
    sumstat_seg_sites()

  stats <- scrm$simulate(model)
  expect_equal(nrow(stats$seg_sites[[1]]), 5)
  expect_equal(nrow(stats$seg_sites[[2]]), 5)
})
