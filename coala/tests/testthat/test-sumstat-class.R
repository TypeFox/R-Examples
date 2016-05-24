context("SumStat Class")

test_that("Sumstat initialization works", {
  ss1 <- sumstat_class$new("1", I)
  expect_true(is.sum_stat(ss1))
  expect_error(ss1$calculate(NULL, NULL, NULL, NULL))

  ss1 <- sumstat_class$new("1", I)
  expect_true(is.sum_stat(ss1))

  expect_false(is.sum_stat(1:3))
})


test_that("getting the name works", {
  ss1 <- sumstat_class$new("test", identity)
  expect_equal(ss1$get_name(), "test")
})


test_that("transformations works", {
  expect_equal(sumstat_class$new("sum", function(x) 1)$transform(1:10), 1)
  expect_equal(sumstat_class$new("sum", sum)$transform(1:3), 6)
  expect_error(sumstat_class$new("sum", 1))
})


test_that("Adding Sumstats to a model works", {
  model <- coal_model(5:6, 10) + sumstat_class$new("1", I)
  expect_equal(get_summary_statistics(model)[[1]]$get_name(), "1")
  expect_error(model + sumstat_class$new("1"))

  model <- model + sumstat_class$new("2", I)
  expect_equal(names(get_summary_statistics(model)), c("1", "2"))
})


test_that("Calculation of sumstats works", {
  stat_sum_class <- R6::R6Class("Stat_Sum", inherit = sumstat_class,
    public = list(calculate = function(seg_sites, trees, files, model) {
      sapply(seg_sites, sum)
    })
  )

  model <- coal_model(5:6, 10) + stat_sum_class$new("sum", identity)
  expect_equal(calc_sumstats(model, list(1:3, 1:5, 7)),
               list(sum = c(6, 15, 7)))
  expect_equal(calc_sumstats(model, list(1:3, 1:5, 7), blub = "test"),
               list(sum = c(6, 15, 7), blub = "test"))
  expect_equal(calc_sumstats(model, list(1:3, 1:5, 7), bla = 1:5),
               list(sum = c(6, 15, 7), bla = 1:5))
})


test_that("Calculation of sumstats from simresults works", {
  stat_sum_class <- R6::R6Class("Stat_Sum", inherit = sumstat_class,
    public = list(calculate = function(seg_sites, trees, files, model) {
      sapply(seg_sites, sum)
    })
  )

  model <- coal_model(5:6, 10) + stat_sum_class$new("sum", identity) +
    stat_sum_class$new("sum2", identity)
  stats <- calc_sumstats_from_sim(list(1:3, 1:5, 7), NULL, "", model, 1:2,
                                  1:3, get_simulator("scrm"))
  expect_equal(stats$sum,  c(6, 15, 7))
  expect_equal(stats$sum2, c(6, 15, 7))
  expect_equal(stats$cmds, 1:3)
  expect_equal(stats$pars, 1:2)
})


test_that("Calculation of sumstats respects transformations", {
  stat_sum_class <- R6::R6Class("Stat_Sum", inherit = sumstat_class,
    public = list(calculate = function(seg_sites, trees, files, model) {
      sapply(seg_sites, sum)
    })
  )

  model <- coal_model(5:6, 10) + stat_sum_class$new("sum", sum)
  expect_equal(calc_sumstats(model, list(1:3, 1:5, 7)), list(sum = 28))

  model <- model + stat_sum_class$new("sum2", function(x) 1:15)
  expect_equal(calc_sumstats(model, list(1:3, 1:5, 7)),
               list(sum = 28, sum2 = 1:15))
})


test_that("summary statistics can be calculated from real data", {
  segsites_list <- list(create_test_segsites(),
                        create_test_segsites(),
                        create_test_segsites())

  model <- coal_model(3, 3, 100) + feat_mutation(5) + sumstat_sfs()

  expect_equal(calc_sumstats_from_data(model, segsites_list),
               list(sfs = c(6, 9)))

  # Test that an error is thrown for missing data
  segsites_list[[2]] <- NA
  expect_error(calc_sumstats_from_data(model, segsites_list))
  segsites_list[[2]] <- create_test_segsites()
  expect_error(calc_sumstats_from_data(coal_model(3, 5), segsites_list))

  # Test import with trios
  model <- coal_model(c(2, 1)) + locus_trio(number = 2) + coala::sumstat_sfs()
  stats <- calc_sumstats_from_data(model, segsites_list, trios = list(1:3, 1:3))
  expect_equal(stats, list(sfs = c(4, 6)))

  model <- coal_model(3) + locus_trio() + locus_single(5) + coala::sumstat_sfs()
  stats <- calc_sumstats_from_data(model, segsites_list, trios = list(1:3, 2))
  expect_equal(stats, list(sfs = c(4, 6)))

  expect_error(calc_sumstats_from_data(model, tree_list = list("12345")))
})
