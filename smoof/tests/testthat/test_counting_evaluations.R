context("counting evaluations")

test_that("countingWrapper counts correctly", {
  fn = makeBBOBFunction(fid = 1L, iid = 1L, dimension = 8L)
  expect_false(doesCountEvaluations(fn))
  fn = addCountingWrapper(fn)
  expect_true(doesCountEvaluations(fn))

  # no evaluations yet
  expect_equal(getNumberOfEvaluations(fn), 0L)

  # now perform 10 function evaluations
  par.set = getParamSet(fn)
  par.mat = matrix(NA, nrow = 8L, ncol = 10L)
  for (i in 1:10) {
    par = unlist(sampleValue(par.set))
    par.mat[, i] = par
    fn(par)
  }

  # we now should have counted 10 function evaluations
  expect_equal(getNumberOfEvaluations(fn), 10L)

  # now perform one call with a matrix of 10 parameter values
  fn(par.mat)

  # and we should have counted 10 more evaluations
  expect_equal(getNumberOfEvaluations(fn), 20L)

  # reset counter works
  resetEvaluationCounter(fn)
  expect_equal(getNumberOfEvaluations(fn), 0L)
})
