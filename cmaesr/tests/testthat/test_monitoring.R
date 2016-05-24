context("CMA-ES monitoring")

test_that("CMA-ES monitoring works well", {
	fn = makeSphereFunction(2L)
  expect_output(cmaes(fn, monitor = makeSimpleMonitor()), "Iteration")
})
