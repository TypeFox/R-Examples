context("Simulator msms")

test_that("calling msms works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  msms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  set.seed(17)
  out_file <- msms$call_msms(msms.args)
  set.seed(17)
  out_file_2 <- msms$call_msms(msms.args)
  set.seed(20)
  out_file_3 <- msms$call_msms(msms.args)
  expect_equal(file.info(out_file_2)$size, file.info(out_file)$size)
  expect_true(file.info(out_file)$size != file.info(out_file_3)$size)
  unlink(c(out_file, out_file_2, out_file_3))
})


test_that("generating msms options works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  model <- coal_model(10, 2) + feat_mutation(5)
  expect_equal(msms$get_cmd(model), "msms 10 2 -t 5 -threads 1 ")
})


test_that("msms can simulate seg. sites", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")

  # Generating Seg. Sites
  model <- coal_model(10, 2, 100) + feat_mutation(5) + sumstat_seg_sites()
  set.seed(110); stats_1 <- msms$simulate(model)
  set.seed(110); stats_2 <- msms$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$seg_sites, is_a("list"))
  expect_equal(length(stats_1$seg_sites), 2)
  expect_equal(length(get_positions(stats_1$seg_sites[[1]])),
               ncol(stats_1$seg_sites[[1]]))

  # With recombination
  model <- model + feat_recombination(1)
  stats <- msms$simulate(model)
  expect_that(stats_1$seg_sites, is_a("list"))
  expect_equal(length(stats_1$seg_sites), 2)
})


test_that("msms can simulate trees", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")

  model <- coal_model(10, 2, 100) + sumstat_trees()
  set.seed(110); stats_1 <- msms$simulate(model)
  set.seed(110); stats_2 <- msms$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
  for (i in 1:2) expect_equal(length(stats_1$trees[[i]]), 1)

  # With inter locus variation
  model <- coal_model(10, 2, 100) +
    feat_selection(strength_A = par_variation(100, 5),
                   population = 1, time = 0.05) +
    sumstat_trees()
  stats <- msms$simulate(model)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
})


test_that("msms can simulate files", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")

  # Generating Files
  test_folder <- tempfile("msms_test_folder")
  model <- coal_model(10, 2, 100) +
    feat_mutation(5) +
    sumstat_file(test_folder)
  stats <- msms$simulate(model)
  expect_true(file.exists(test_folder))
  expect_true(file.exists(stats$file[1]))
  unlink(test_folder, recursive = TRUE)
})


test_that("msms_simulate works with inter-locus variation", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")

  model_tmp <- coal_model(5, 2) +
    feat_mutation(par_variation(par_range("theta", 1, 5), 17)) +
    sumstat_seg_sites()
  expect_true(has_variation(model_tmp))

  sum_stats <- msms$simulate(model_tmp, c(theta = 3))
  expect_is(sum_stats$seg_sites, "list")
  expect_equal(length(sum_stats$seg_sites), 2)
})


test_that("simulating unphased data works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  model <- coal_model(5, 2, ploidy = 2) +
    feat_unphased(1) +
    feat_mutation(4, fixed_number = TRUE) +
    sumstat_seg_sites()

  stats <- msms$simulate(model)
  expect_equal(nrow(stats$seg_sites[[1]]), 5)
})


test_that("msms can simulate locus trios", {
  if (!has_msms()) skip("msms not installed")
  stat <- get_simulator("msms")$simulate(model_trios())
  expect_that(get_trio_locus(stat$seg_sites[[1]]), is_a("numeric"))
  expect_true(all(get_trio_locus(stat$seg_sites[[1]]) %in% -1:1))
  expect_true(all(get_positions(stat$seg_sites[[1]]) >= 0))
  expect_true(all(get_positions(stat$seg_sites[[1]]) <= 1))
})


test_that("msms can be added manually", {
  if (!has_msms()) skip("msms not installed")
  msms_jar <- get_simulator("msms")$get_info()["jar"]
  java <- get_simulator("msms")$get_info()["java"]
  activate_msms(msms_jar, java, 199)
  expect_equal(get_simulator("msms")$get_priority(), 199)
  expect_error(use_msms(tempfile("not-existant"), tempfile("not-existant")))
})


test_that("msms supports size changes in one pop models", {
  if (!has_msms()) skip("msms not installed")

  model <- coal_model(40, 1) +
    feat_mutation(1) +
    feat_size_change(0.1, population = 1, time = 0.1) +
    sumstat_sfs()

  stat <- get_simulator("msms")$simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})


test_that("msms supports growth in one pop models", {
  if (!has_msms()) skip("msms not installed")

  model <- coal_model(40, 1) +
    feat_mutation(1) +
    feat_growth(0.1, population = 1, time = 0.1) +
    sumstat_sfs()

  stat <- get_simulator("msms")$simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})
