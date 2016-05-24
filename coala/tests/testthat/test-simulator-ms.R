context("Simulator ms")


test_that("parsing positions works", {
  expect_equal(parse_ms_positions("positions: 0.0010 0.0474 0.3171"),
               c(0.001, 0.0474, 0.3171))
  expect_equal(parse_ms_positions("positions: 0.1 0.2 0.3 0.4 0.5"), 1:5 / 10)
  expect_equal(parse_ms_positions("positions: 0.1"), 0.1)

  expect_error(parse_ms_positions("0.1 0.2 0.3"))
  expect_error(parse_ms_positions(" "))
  expect_error(parse_ms_positions("segsites: 0"))
})


test_that("Parsing ms output works", {
  sim_output <- tempfile("sim_output")

  cat("ms 3 1 -t 1 -r 1 20 -T
30461 15911 34727

//
[2](2:0.865,(1:0.015,3:0.015):0.850);
[3](2:0.865,(1:0.015,3:0.015):0.850);
[4](2:1.261,(1:0.015,3:0.015):1.246);
[11](2:1.261,(1:0.015,3:0.015):1.246);
segsites: 5
positions: 0.2046 0.2234 0.2904 0.6209 0.9527
01100
10011
01100

//
[2](3:0.613,(1:0.076,2:0.076):0.537);
[18](3:0.460,(1:0.076,2:0.076):0.384);
segsites: 2
positions: 0.3718 0.8443
01
01
10
", file = sim_output);

  output <- parse_ms_output(list(sim_output), 3, 2)


  ss1 <- create_segsites(matrix(c(0, 1, 1, 0, 0,
                                  1, 0, 0, 1, 1,
                                  0, 1, 1, 0, 0), 3, 5, TRUE),
                         c(0.2046, 0.2234, 0.2904, 0.6209, 0.9527))

  ss2 <- create_segsites(matrix(c(0, 0, 1, 1, 1, 0), 3, 2), c(0.3718, 0.8443))

  trees1 <- c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
              "[3](2:0.865,(1:0.015,3:0.015):0.850);",
              "[4](2:1.261,(1:0.015,3:0.015):1.246);",
              "[11](2:1.261,(1:0.015,3:0.015):1.246);")
  trees2 <- c("[2](3:0.613,(1:0.076,2:0.076):0.537);",
              "[18](3:0.460,(1:0.076,2:0.076):0.384);")

  expect_equal(output, list(segsites = list(ss1, ss2),
                            trees = list(trees1, trees2)))

  expect_error(parse_ms_output(list(sim_output), 3, 1))
  expect_error(parse_ms_output(list(sim_output), 3, 3))

  output <- parse_ms_output(list(sim_output, sim_output), 3, 4)
  expect_equal(output, list(segsites = list(ss1, ss2, ss1, ss2),
                            trees = list(trees1, trees2, trees1, trees2)))

  output2 <- parse_ms_output(list(c(sim_output, sim_output)), 3, 4)
  expect_equal(output, output2)

  output <- parse_ms_output(list(c(sim_output, sim_output),
                                 c(sim_output, sim_output)), 3, 8)
  expect_equal(length(output$segsites), 8)
  expect_equal(length(output$trees), 8)

  expect_error(parse_ms_output(tempfile("test_ms_out"), 4, 4))

  unlink(sim_output)
})


test_that("the ms sim program exists", {
  if (!has_ms()) skip("ms not installed")
  expect_that(get_simulator("ms"), is_a("simulator"))
})


test_that("msms can simulate seg. sites", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  # Generating Seg. Sites
  model <- coal_model(10, 2, 100) + feat_mutation(5) + sumstat_seg_sites()
  set.seed(110); stats_1 <- ms$simulate(model)
  set.seed(110); stats_2 <- ms$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$seg_sites, is_a("list"))
  expect_equal(length(stats_1$seg_sites), 2)
  expect_equal(length(get_positions(stats_1$seg_sites[[1]])),
               ncol(stats_1$seg_sites[[1]]))

  # With recombination
  model <- model + feat_recombination(1)
  stats <- ms$simulate(model)
  expect_that(stats_1$seg_sites, is_a("list"))
  expect_equal(length(stats_1$seg_sites), 2)
})


test_that("ms can simulate trees", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  model <- coal_model(10, 2, 100) + sumstat_trees()
  set.seed(110); stats_1 <- ms$simulate(model)
  set.seed(110); stats_2 <- ms$simulate(model)
  expect_equal(stats_1, stats_2)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
  for (i in 1:2) expect_equal(length(stats_1$trees[[i]]), 1)

  # With inter locus variation
  model <- coal_model(10, 2, 100) +
    feat_recombination(par_variation(5, 2)) +
    sumstat_trees()
  stats <- ms$simulate(model)
  expect_that(stats_1$trees, is_a("list"))
  expect_equal(length(stats_1$trees), 2)
})


test_that("msms can simulate files", {
  if (!has_ms()) skip("msms not installed")
  ms <- get_simulator("ms")

  # Generating Files
  test_folder <- tempfile("ms_test_folder")
  model <- coal_model(10, 2, 100) +
    feat_mutation(5) +
    sumstat_file(test_folder)
  stats <- ms$simulate(model)
  expect_true(file.exists(test_folder))
  expect_true(file.exists(stats$file[1]))
  unlink(test_folder, recursive = TRUE)
})


test_that("Saving the simulation cmds works", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  model <- model_theta_tau() + locus_single(15)
  stats <- ms$simulate(model, c(tau = 1, theta = 10))
  expect_that(stats$cmds, is_a("character"))
  expect_equal(length(stats$cmds), 2)
  expect_true(grepl("^ms 25 10 ", stats$cmds[1]))
  expect_equal(length(stats$cmds[[2]]), 1)
  expect_true(grepl("^ms 25 1 ", stats$cmds[2]))

  model <- coal_model(5, 10) +
    locus_single(15) +
    feat_mutation(1) +
    feat_recombination(par_zero_inflation(5, .5)) +
    sumstat_sfs()
  stats <- ms$simulate(model)
  expect_equal(length(stats$cmds), 3)
})


test_that("simulating unphased data works", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  model <- coal_model(3, 2, ploidy = 2) +
    feat_unphased(1) +
    feat_mutation(4, fixed_number = TRUE) +
    sumstat_seg_sites()

  stats <- ms$simulate(model)
  expect_equal(nrow(stats$seg_sites[[1]]), 3)
  expect_equal(nrow(stats$seg_sites[[2]]), 3)
})


test_that("ms can simulate locus trios", {
  if (!has_ms()) skip("ms not installed")
  stat <- get_simulator("ms")$simulate(model_trios())

  expect_that(get_trio_locus(stat$seg_sites[[1]]), is_a("numeric"))
  expect_true(all(get_trio_locus(stat$seg_sites[[1]]) %in% -1:1))
  expect_true(all(get_positions(stat$seg_sites[[1]]) >= 0))
  expect_true(all(get_positions(stat$seg_sites[[1]]) <= 1))
})


test_that("ms works with scientific notation", {
  if (!has_ms()) skip("ms not installed")
  ms <- get_simulator("ms")

  model <- coal_model(5, 1, 1e8) + feat_recombination(1)
  template <- ms$create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))

  model <- coal_model(5, 1, 1000) + feat_recombination(1e8)
  template <- ms$create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))

  model <- coal_model(5, 1, 1000) + feat_mutation(1e8)
  template <- ms$create_cmd_tempalte(model)
  opts <- fill_cmd_template(template, model, numeric(0), 1)
  expect_true(grepl("100000000", opts$command[1]))
})


test_that("ms can simulate zero inflation", {
  if (!has_ms()) skip("ms not installed")

  model <- model_theta_tau() +
    feat_recombination(par_zero_inflation(1, .5)) +
    locus_averaged(4, 100) +
    locus_single(10)
  stats <- get_simulator("ms")$simulate(model, c(tau = 1, theta = 5))
  expect_that(stats, is_a("list"))
})
