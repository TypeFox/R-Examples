context("Model Class")


test_that("creating models works", {
  model <- coal_model(11:12, 111, 1234)
  expect_true(is.model(model))
  expect_equal(get_sample_size(model), 11:12)
  expect_equal(get_locus_number(model), 111)
  expect_equal(get_locus_length(model, 1), 1234)

  expect_is(model$id, "character")
})


test_that("adding parameters works", {
  model <- coal_model(11:12, 100) + par_range("p1", 1, 5)
  par_table <- get_parameter_table(model)
  expect_equal("p1", par_table$name)
  expect_equal(1, par_table$lower.range)
  expect_equal(5, par_table$upper.range)

  expect_that(get_parameter(model), is_a("list"))
  expect_equal(length(get_parameter(model)), 1)
  expect_true(is.par(get_parameter(model)[[1]]))

  model <- model + par_range("p2", 1, 5)
  expect_equal(length(get_parameter(model)), 2)

  test <- list(1:10)
  class(test) <- "BLUB"
  expect_error(model + test)
})


test_that("adding features works", {
  expect_equal(length(get_features(coal_model(5))), 1)
  model <- coal_model(5) + feature_class$new(1, 1, 5)
  expect_equal(length(get_features(model)), 2)
  model <- model + feature_class$new(2, 1, 3)
  expect_equal(length(get_features(model)), 3)
})


test_that("test get_summary_statistics", {
  expect_equal(get_summary_statistics(coal_model(1:2)), list())
  expect_equal(length(get_summary_statistics(model_theta_tau())),  1)
  expect_true(is.sum_stat(get_summary_statistics(model_theta_tau())[[1]]))
})


test_that("test that scaling of model works", {
  model <- coal_model(11:12, 10) +
    locus_averaged(24, 10) +
    locus_averaged(25, 15) +
    locus_single(101) +
    locus_single(102)
  model <- scale_model(model, 5)

  expect_equal(get_locus_number(model), 14)
  expect_equal(get_locus_number(model, 1), 2)
  expect_equal(get_locus_number(model, 2), 5)
  expect_equal(get_locus_number(model, 5), 1)
})


test_that("get loci length and number works", {
  model <- coal_model(10, 11, 101) +
    locus_averaged(12, 102) +
    locus_trio(locus_length = 1:3, distance = 10:11)

  expect_equal(get_locus_number(model), 24)
  expect_equal(get_locus_length(model, 1), 101)
  expect_equal(get_locus_length(model, 5), 101)
  expect_equal(get_locus_length(model, 11), 101)
  expect_equal(get_locus_length(model, 15), 102)
  expect_equal(get_locus_length(model, 23), 102)
  expect_equal(get_locus_length(model, 24), 27)

  expect_equal(get_locus_length(model, group = 1), 101)
  expect_equal(get_locus_length(model, group = 2), 102)
  expect_equal(get_locus_length(model, group = 3), 27)

  expect_equivalent(get_locus_length(model, 1, total = FALSE), 101)
  expect_equivalent(get_locus_length(model, 24, total = FALSE),
                    c(1, 10, 2, 11, 3))

  expect_equal(get_locus_length(model, total = TRUE), c(101, 102, 27))
  expect_equivalent(get_locus_length(model, total = FALSE),
                    matrix(c(0, 0, 101, 0, 0,
                             0, 0, 102, 0, 0,
                             1, 10, 2, 11, 3), 3, byrow = TRUE))
})


test_that("locus length matrix generations works", {
  # Multiple loci with equal length
  expect_equivalent(get_locus_length_matrix(model_theta_tau()),
                    matrix(c(0, 0, 1000, 0, 0, 10), 1, 6))

  # Multiple loci with differnt length
  model <- model_theta_tau() +
    locus_single(21) +
    locus_single(22) +
    locus_single(23)

  expect_equivalent(get_locus_length_matrix(model),
                    matrix(c(0, 0, 1000, 0, 0, 10,
                             0, 0, 21, 0, 0, 1,
                             0, 0, 22, 0, 0, 1,
                             0, 0, 23, 0, 0, 1), 4, 6, TRUE))

  # Test with scaling
  model <- scale_model(model, 5)
  expect_equivalent(get_locus_length_matrix(model),
                    matrix(c(0, 0, 1000, 0, 0, 2,
                             0, 0, 21, 0, 0, 1,
                             0, 0, 22, 0, 0, 1,
                             0, 0, 23, 0, 0, 1), 4, 6, TRUE))
})


test_that("getting the available Populations works", {
  model <- coal_model(10:11, 100)
  expect_equal(get_populations(model), 1:2)
  expect_equal(get_populations(model_theta_tau()), 1:2)
  expect_equal(get_populations(model_hky()), 1:3)
  expect_equal(get_populations(model + feat_sample(1:5)), 1:5)
})


test_that("get population individuals works", {
  expect_equal(get_population_individuals(model_theta_tau(), 1), 1:10)
  expect_equal(get_population_individuals(model_theta_tau(), 2), 11:25)
  expect_equal(get_population_individuals(model_theta_tau(), "all"), 1:25)
  expect_error(get_population_individuals(model_theta_tau(), 3))
  expect_error(get_population_individuals(model_theta_tau(), "al"))

  # With an outgroup
  expect_equal(get_population_individuals(model_hky(), "all"), 1:6)
  expect_equal(get_population_individuals(model_hky(), 1), 1:3)
  expect_equal(get_population_individuals(model_hky(), 2), 4:6)
  expect_error(get_population_individuals(model_hky(), 3))

  model <- coal_model(1:3) + feat_outgroup(2)
  expect_equal(get_population_individuals(model, "all"), 1:4)
  expect_equal(get_population_individuals(model, 1), 1)
  expect_error(get_population_individuals(model, 2))
  expect_equal(get_population_individuals(model, 3), 2:4)
})


test_that("getting indiviuals in polyploid models works", {
  model <- coal_model(1:5, ploidy = 2)
  expect_equal(get_population_individuals(model, 1, haploids = FALSE), 1)
  expect_equal(get_population_individuals(model, 2, haploids = FALSE), 2:3)
  expect_equal(get_population_individuals(model, 3, haploids = FALSE), 4:6)
  expect_equal(get_population_individuals(model, 4, haploids = FALSE), 7:10)
  expect_equal(get_population_individuals(model, 5, haploids = FALSE), 11:15)
})


test_that("getting the ploidy and individuals works", {
  model <- model_theta_tau()
  expect_equal(get_ploidy(model), 1L)
  expect_equal(get_samples_per_ind(model), 1L)
  sample_size <- get_sample_size(model)
  expect_false(is_unphased(model))

  model <- coal_model(sample_size, ploidy = 4) + feat_unphased(2)
  expect_equal(get_ploidy(model), 4L)
  expect_equal(get_samples_per_ind(model), 2L)
  expect_equal(get_sample_size(model), sample_size * 2)
  expect_equal(get_sample_size(model, for_sim = TRUE), sample_size * 4)
  expect_true(is_unphased(model))
})


test_that("print works on models", {
  # Printing an empty model works
  out <- capture.output(print(coal_model(5)))
  expect_that(length(out), is_more_than(0))

  # Printing parameters works
  out <- capture.output(print(coal_model(5) + par_range("abc", 1, 5)))
  expect_that(length(grep("abc", out)), is_more_than(0))

  # Printing loci works
  out <- capture.output(print(coal_model(5) + locus_single(3131)))
  expect_that(length(grep("3131", out)), is_more_than(0))

  out <- capture.output(print(coal_model(5) +
                                locus_single(3131) +
                                locus_single(3131)))
  expect_that(length(grep("3131", out)), is_more_than(0))
})


test_that("getting par names works", {
  expect_equal(get_par_names(coal_model(5)), character(0))

  model <- coal_model(5) + par_range("a", 1, 2) + par_range("b", 2, 3)
  expect_equal(get_par_names(model), c("a", "b"))
  expect_equal(get_par_names(model, TRUE), c("a", "b"))

  model <- model + par_prior("c", 1)
  expect_equal(get_par_names(model), c("a", "b", "c"))
  expect_equal(get_par_names(model, TRUE), c("a", "b"))
})


test_that("getting model command works", {
  cmd <- get_cmd(model_theta_tau())
  expect_that(cmd, is_a("character"))
  expect_that(nchar(cmd), is_more_than(0))
})


test_that("has_trios works", {
  expect_false(has_trios(model_theta_tau()))
  expect_false(has_trios(model_gtr()))
  expect_true(has_trios(model_trios()))
})


test_that("creating a parameter table works ", {
  expect_equal(get_parameter_table(coal_model(5)),
               data.frame(name = character(0),
                          lower.range = numeric(0),
                          upper.range = numeric(0),
                          stringsAsFactors = FALSE))

  model <- coal_model(5:6, 10, 100) + par_range("theta", 1, 2)
  expect_equal(get_parameter_table(model),
               data.frame(name = "theta", lower.range = 1, upper.range = 2,
                          stringsAsFactors = FALSE))

  model <- coal_model(5:6, 10, 100) +
    par_range("theta", 1, 2) +
    par_range("tau", 5, 6)
  expect_equal(get_parameter_table(model),
               data.frame(name = c("theta","tau"),
                          lower.range = c(1, 5),
                          upper.range = c(2, 6),
                          stringsAsFactors = FALSE))

  expect_error(get_parameter_table(model + par_prior("x", rnorm(1))))
})


test_that("model checking give not errors", {
  capture.output(check_model(model_theta_tau()))
  capture.output(check_model(model_gtr()))
  capture.output(check_model(model_hky()))
  capture.output(check_model(model_trios()))
})


test_that("model parts can be combined into a partial model", {
  incomplete_model <- feat_growth(1, 1) + feat_growth(2, 2)
  expect_true(is_partial_model(incomplete_model))
  expect_equal(length(incomplete_model), 2)

  incomplete_model <- sumstat_sfs() + sumstat_dna()
  expect_true(is_partial_model(incomplete_model))
  expect_equal(length(incomplete_model), 2)

  incomplete_model <- par_const(5) + par_const(7)
  expect_true(is_partial_model(incomplete_model))
  expect_equal(length(incomplete_model), 2)

  incomplete_model <- locus_single(1) + locus_averaged(2, 10)
  expect_true(is_partial_model(incomplete_model))
  expect_equal(length(incomplete_model), 2)
})


test_that("partial models can be extended", {
  incomplete_model <- feat_growth(1, 1) +
    feat_mutation(5) +
    feat_recombination(7)
  expect_true(is_partial_model(incomplete_model))
  expect_equal(length(incomplete_model), 3)
  expect_equal(incomplete_model[[1]], feat_growth(1, 1))
  expect_equal(incomplete_model[[2]], feat_mutation(5))
  expect_equal(incomplete_model[[3]], feat_recombination(7))

  incomplete_model_2 <- incomplete_model + sumstat_sfs()
  expect_true(is_partial_model(incomplete_model_2))
  expect_equal(length(incomplete_model_2), 4)
  expect_equal(incomplete_model_2[[4]], sumstat_sfs())
})


test_that("partial models can be added to models", {
  incomplete_model <- feat_growth(1, 1) +
    feat_mutation(5) +
    feat_recombination(par_const(8)) +
    locus_averaged(10, 100) +
    par_const(6) +
    sumstat_sfs()
  model <- coal_model(10) + incomplete_model

  model_direct <- coal_model(10) +
    feat_growth(1, 1) +
    feat_mutation(5) +
    feat_recombination(par_const(8)) +
    locus_averaged(10, 100) +
    par_const(6) +
    sumstat_sfs()
  model_direct$id <- model$id

  expect_equal(model, model_direct)
})


test_that("printing partical models works", {
  incomplete_model <- feat_growth(1, 1) +
    feat_mutation(5) +
    feat_recombination(par_const(8))
  expect_output(print(incomplete_model), "growth")
  expect_output(print(incomplete_model), "Mutation")
  expect_output(print(incomplete_model), "Recombination")
})
