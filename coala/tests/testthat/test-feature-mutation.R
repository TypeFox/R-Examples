context("Feature Mutation")


test_that("Creation of finite sites features works", {
  expect_error(feat_mutation(5, model = "BLUB"))
  expect_error(feat_mutation(5, model = "F84"))

  # HKY
  feat <- feat_mutation(17, "HKY", tstv_ratio = 2,
                        base_frequencies = rep(.25, 4))
  expect_equal(feat$get_tstv_ratio(), 2)
  expect_equal(feat$get_base_frequencies(), rep(.25, 4))
  expect_error(feat_mutation(5, model = "HKY"))
  expect_error(feat_mutation(5, model = "HKY", tstv_ratio = 2))
  expect_error(feat_mutation(5, model = "HKY", base_frequencies = rep(.25, 4)))
  expect_error(feat_mutation(5, model = "HKY", tstv_ratio = 2,
                             base_frequencies = rep(.5, 4)))


  # GTR rates
  feat <- feat_mutation(5, model = "GTR", gtr_rates = 1:6)
  expect_equal(feat$get_gtr_rates(), 1:6)
  expect_error(feat_mutation(5, model = "GTR", gtr_rates = 1:5))
  expect_error(feat_mutation(5, model = "GTR", gtr_rates = "Blub"))
})


test_that("Parsing mutation to scrm args works", {
  feat <- feat_mutation(5)
  scrm_arg <- conv_to_scrm_arg(feat, NULL)
  expect_that(scrm_arg, is_a("character"))
  expect_true(grepl("-t", scrm_arg))
  expect_true(grepl("5", scrm_arg))

  scrm <- get_simulator("scrm")
  model <- coal_model(15, 1) + feat_mutation(par_range("theta", 1, 2))
  expect_equal(scrm$get_cmd(model), "scrm 15 1 -t theta ")
  model <- coal_model(15, 1) + feat_mutation(5)
  expect_equal(scrm$get_cmd(model), "scrm 15 1 -t 5 ")
  model <- coal_model(15, 1) +
    par_range("x", 1, 2) +
    feat_mutation(par_expr(2 * x))
  expect_equal(scrm$get_cmd(model), "scrm 15 1 -t 2 * x ")

  expect_error(conv_to_scrm_arg(feat_mutation(5, "GTR"), NULL))
})


test_that("Parsing mutation to seqgen args works", {
  if (!has_seqgen()) skip("seqgen not installed")
  sg <- get_simulator("seqgen")
  model <- coal_model(10, 1) + feat_mutation(5, model = "GTR", gtr_rates = 1:6)
  sg$get_cmd(model)
})


test_that("using a fixed number of mutations works with ms", {
  if (!has_ms()) skip("ms not installed")
  model <- coal_model(5, 1) +
    feat_mutation(5, fixed_number = TRUE) +
    sumstat_sfs()
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 5 1 -s 5 ")
})


test_that("using a fixed number of mutations works with msms", {
  if (!has_msms()) skip("msms not installed")
  model <- coal_model(5, 1) +
    feat_mutation(5, fixed_number = TRUE) +
    sumstat_sfs()
  expect_true(grepl("-s 5", get_simulator("msms")$get_cmd(model)))
})


test_that("scrm rejects models with a fixed number of mutations", {
  expect_error(get_simulator("scrm")$get_cmd(model))
})


test_that("seq-gen rejects models with a fixed number of mutations", {
  if (!has_seqgen()) skip("seqgen not installed")
  model <- coal_model(5, 1) +
    feat_mutation(5, fixed_number = TRUE, model = "GTR", gtr_rates = 1:6) +
    sumstat_sfs()
  expect_error(get_simulator("seqgen")$get_cmd(model))
})
