context("CMA-ES stopping conditions")

test_that("CMA-ES stopping conditions work fine", {
  # test everything on sphere function here
  fn = makeSphereFunction(2L)

  # stop on maximum iterations reached
  max.iters = 50L
  control = list(stop.ons = list(stopOnMaxIters(max.iters)))
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(grepl("iterations", res$message, ignore.case = TRUE))
  expect_equal(res$n.iters, max.iters)

  # test maximum time
  max.time = 2
  control = list(stop.ons = list(stopOnTimeBudget(max.time)))
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(grepl("budget", res$message, ignore.case = TRUE))
  expect_true((res$past.time - max.time) < 1)

  # stop on low gap to optimal params
  opt.param = getGlobalOptimum(fn)$param
  control = list(stop.ons = list(stopOnOptParam(as.numeric(opt.param))))
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(grepl("parameters", res$message, ignore.case = TRUE))

  # stop on maximum number of function evaluations
  max.evals = 100L
  control = list(stop.ons = list(stopOnMaxEvals(max.evals)))
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(grepl("evaluations", res$message, ignore.case = TRUE))

  # stop on low gap to optimal objective function value
  opt.value = getGlobalOptimum(fn)$value
  control = list(stop.ons = list(stopOnOptValue(opt.value)))
  res = cmaes(fn, control = control, monitor = NULL)
  expect_true(grepl("function value", res$message, ignore.case = TRUE))
})
