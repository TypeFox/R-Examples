context("getIndex")

test_that("getIndex", {
  r = makeTestRegistry()
  p1 = addProblem(r, "one", 1)
  p2 = addProblem(r, "two", 2)
  a1 = addAlgorithm(r, "A", fun=function(static, dynamic) 1)
  a2 = addAlgorithm(r, "B", fun=function(static, dynamic) 1)
  addExperiments(r, list(makeDesign(p1), makeDesign(p2)), list(makeDesign(a1), makeDesign(a2)), repls=2)

  expect_equal(getIndex(r, by.prob=TRUE),
               list(prob = factor(rep(c("one", "two"), each = 4))))
  expect_equal(getIndex(r, by.algo=TRUE),
               list(algo = factor(rep(rep(c("A", "B"), each = 2), 2))))
  expect_true(length(getIndex(r, by.algo=TRUE, by.prob=TRUE)) == 2L)
  expect_true(length(getIndex(r, ids=integer(0))) == 0)
  expect_true(length(getIndex(r, ids=integer(0), by.prob=TRUE)) == 1)

  r = makeTestRegistry()
  p1 = addProblem(r, "one", dynamic = function(i, k) i^k)
  a1 = addAlgorithm(r, "A", fun=function(dynamic, k) dynamic^k)
  addExperiments(r,
                 makeDesign("one", exhaustive = list(i =1:3, k = 2)),
                 makeDesign("A", exhaustive = list(k=1:3)))
  foo = 1
  expect_equal(getIndex(r, by.prob.pars=i == foo)[[1]],
               factor(c(rep(TRUE, 3), rep(FALSE, 6)), levels=c("FALSE", "TRUE")))
  expect_true(grepl("i == k", names(getIndex(r, by.prob.pars = i == k))))
  expect_true(length(getIndex(r, by.prob.pars = k, by.algo.pars = k)) == 2)
})
