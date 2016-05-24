context("Interface to abc")

test_that("it creates the abc parameter data.frame", {
  sim_results <- simulate(model_theta_tau(), pars = c(1, 5), nsim = 3)
  sim_results[[2]]$pars <- c(tau = 2, theta = 10)
  sim_results[[3]]$pars <- c(tau = 3, theta = 15)

  expect_equal(create_abc_param(sim_results),
               data.frame(tau = 1:3, theta = 1:3 * 5))

  sim_results <- simulate(model_theta_tau(), pars = c(1, 5), nsim = 1)
  expect_equal(create_abc_param(sim_results), data.frame(tau = 1, theta = 5))
})


test_that("it creates abc sumstat data.frame", {
  model <- model_theta_tau()
  sim_results <- simulate(model, pars = c(1, 5), nsim = 2)
  expect_warning(abc_sumstat <- create_abc_sumstat(sim_results, model))
  expect_equivalent(abc_sumstat[1, ], as.vector(sim_results[[1]]$jsfs))
  expect_equivalent(abc_sumstat[2, ], as.vector(sim_results[[2]]$jsfs))

  sim_results <- simulate(model, pars = c(1, 5), nsim = 1)
  expect_warning(abc_sumstat <- create_abc_sumstat(sim_results, model))
  expect_equivalent(abc_sumstat[1, ], as.vector(sim_results$jsfs))

  model <- coal_model(10, 1) + feat_mutation(5) +
    sumstat_sfs() + sumstat_mcmf()
  sim_results <- simulate(model, 5)
  expect_warning(abc_sumstat <- create_abc_sumstat(sim_results, model), NA)
  expect_equal(dim(abc_sumstat), c(5, 10))
})
