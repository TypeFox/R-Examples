context("Model Simulation")

test_that("basic models can be simulated", {
  expect_error(simulate(model_theta_tau(), pars = 1))
  expect_error(simulate(model_theta_tau(), pars = 1:3))
  expect_error(simulate(model_theta_tau(), pars = c(2, 50)))
  sum.stats <- simulate(model_theta_tau(), pars = c(1, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)
})


test_that("parallel simulations are reproducible", {
  skip_on_os("windows")
  model <- coal_model(10) + feat_mutation(5) + sumstat_seg_sites() +
    locus_single(10) + locus_single(10) + locus_single(10) + locus_single(10)

  res <- simulate(model, seed = 215, cores = 2)
  expect_false(identical(res[[1]], res[[2]]))

  res2 <- simulate(model, seed = 215, cores = 2)
  expect_equal(res, res2)

  res3 <- simulate(model, seed = 215, cores = 1)
  expect_equal(res, res3)
})


test_that("models with multiple loci groups can be simulated", {
  model <- model_theta_tau() + locus_averaged(2, 100) + sumstat_seg_sites()
  sum_stats <- simulate(model, pars = c(1, 5))
  expect_equal(length(sum_stats$seg_sites), get_locus_number(model))
})


test_that("model without ranged parameters can be simualted", {
  model <- coal_model(5, 1) + feat_mutation(5) + sumstat_sfs("sfs")
  stat <- simulate(model)
  expect_false(is.null(stat$sfs))
})


test_that("models with priors can be simulated", {
  model <- coal_model(5, 1) +
    feat_mutation(5) +
    feat_recombination(par_prior("r", stats::rbinom(1, 3, .5))) +
    sumstat_sfs("sfs")
  stats <- simulate(model)
  expect_equal(names(stats$pars), "r")
  expect_true(all(stats$pars %in% 0:3))
})


test_that("simulating with more than one repetition works", {
  sim <- simulate(model_theta_tau(), 2, seed = 17, pars = c(1, 5))
  expect_that(sim, is_a("list"))
  expect_equal(length(sim), 2)
  expect_equal(sim[[1]]$pars, c(tau = 1, theta = 5))
  expect_equal(sim[[2]]$pars, c(tau = 1, theta = 5))
})
