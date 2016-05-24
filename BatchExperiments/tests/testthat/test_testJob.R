context("testJob")
test_that("testJob", {
  reg = makeTestRegistry()
  p1 = addProblem(reg, "p1", 1)
  a1 = addAlgorithm(reg, id="a1", fun=function(static, dynamic) static)
  addExperiments(reg, p1, a1)
  expect_equal(testJob(reg, 1L, external=FALSE), 1)
  #FIXME: for some reason this test does not run in "make check"
  # is should also be tested as well in BJ, not only here,
  # with external = TRUE!
  if (isExpensiveExampleOk() && interactive())
    expect_equal(testJob(reg, 1L, external=TRUE), 1)
})
