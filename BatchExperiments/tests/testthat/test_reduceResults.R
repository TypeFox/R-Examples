context("reduceResults")

test_that("reduceResults", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addProblem(reg, "p2", static = 2)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) static*1)
  addAlgorithm(reg, id = "a2", fun = function(static, dynamic) static*2)
  addExperiments(reg, c("p1", "p2"), c("a1", "a2"), 2)
  submitJobs(reg)
  waitForJobs(reg)
  data = reduceResults(reg, fun = function(aggr, job, res)
    rbind(aggr, data.frame(id = job$id, prob = job$prob.id, algo = job$algo.id,
      repl = job$repl, y = res, stringsAsFactors = FALSE)),
    init = data.frame()
  )
  data2 = data.frame(
    id = 1:8,
    prob = c("p1", "p1", "p1", "p1", "p2", "p2", "p2", "p2"),
    algo = c("a1", "a1", "a2", "a2", "a1", "a1", "a2", "a2"),
    repl = as.integer(c(1,2,1,2,1,2,1,2)),
    y = c(1,1,2,2,2,2,4,4),
    stringsAsFactors = FALSE
  )
  expect_equal(data, data2)
  data = reduceResultsExperiments(reg, fun = function(job, res) data.frame(y = res), strings.as.factors = FALSE)
  attr(data2, "prob.pars.names") = character(0)
  attr(data2, "algo.pars.names") = character(0)
  class(data2) = c("ReducedResultsExperiments", "data.frame")
  expect_equal(data, data2)

  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic, a) 1*static)
  addAlgorithm(reg, id = "a2", fun = function(static, dynamic, a, b) 2*static)
  ad1 = makeDesign("a1", exhaustive = list(a = 1:2))
  ad2 = makeDesign("a2", exhaustive = list(a = 3, b = c("b1", "b2")))
  addExperiments(reg, "p1", list(ad1, ad2))
  submitJobs(reg)
  waitForJobs(reg)
  data = reduceResultsExperiments(reg, fun = function(job, res) data.frame(y = res), strings.as.factors = TRUE)
  data2 = data.frame(
    id = 1:4,
    prob = c("p1", "p1", "p1", "p1"),
    algo = c("a1", "a1", "a2", "a2"),
    a = c(1,2,3,3),
    repl = rep(1L, 4),
    y = c(1,1,2,2),
    b = c(NA,NA,"b1","b2"),
    stringsAsFactors = TRUE)
  attr(data2, "prob.pars.names") = character(0)
  attr(data2, "algo.pars.names") = c("a", "b")
  class(data2) = c("ReducedResultsExperiments", "data.frame")
  expect_equal(getResultVars(data, "prob.pars"), character(0))
  expect_equal(getResultVars(data, "algo.pars"), c("a", "b"))
  expect_equal(getResultVars(data, "group"), c("prob", "algo", "a", "b"))
  expect_equal(getResultVars(data, "result"), "y")

  reg = makeTestRegistry()
  reg$seed = 1
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic, a) 1*static)
  addExperiments(reg, "p1", ad1)
  submitJobs(reg)
  waitForJobs(reg)
  data = reduceResultsExperiments(reg, fun = function(job, res) data.frame(y = res, seed = job$seed),
    strings.as.factors = FALSE)
  data2 = data.frame(
    id = 1:2,
    prob = c("p1", "p1"),
    algo = c("a1", "a1"),
    a = c(1,2),
    repl = rep(1L, 2),
    y = c(1,1),
    seed = 1:2,
    stringsAsFactors = FALSE)
  attr(data2, "prob.pars.names") = character(0)
  attr(data2, "algo.pars.names") = c("a")
  class(data2) = c("ReducedResultsExperiments", "data.frame")
  expect_equal(data, data2)
})

test_that("reduceResultsExperiments works on empty id sets", {
  reg = makeTestRegistry()
  reg$seed = 1
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic, a) 1*static)
  addExperiments(reg, repls = 10)
  data2 = data.frame()
  attr(data2, "prob.pars.names") = character(0)
  attr(data2, "algo.pars.names") = character(0)
  class(data2) = c("ReducedResultsExperiments", "data.frame")
  expect_equal(
    reduceResultsExperiments(reg, fun = function(job, res) data.frame(y = res, seed = job$seed),
      strings.as.factors = FALSE),
    data2)
  submitJobs(reg)
  waitForJobs(reg)
  expect_equal(
    reduceResultsExperiments(reg, fun = function(job, res) data.frame(y = res, seed = job$seed),
      strings.as.factors = FALSE, ids = integer(0)),
    data2)
})

test_that("params are available in reduceResults", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", dynamic = function(static, alpha) alpha)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic, beta) dynamic*beta)
  pd1 = makeDesign("p1", exhaustive = list(alpha = 1:2))
  ad1 = makeDesign("a1", exhaustive = list(beta = 3:4))
  addExperiments(reg, pd1, ad1)
  submitJobs(reg)
  waitForJobs(reg)
  data = reduceResults(reg, fun = function(aggr, job, res) {
    rbind(aggr, data.frame(alpha = job$prob.pars$alpha, beta = job$algo.pars$beta, y = res))
  }, init = data.frame())
  expect_equal(data, data.frame(
    alpha = c(1, 1, 2, 2),
    beta = c(3, 4, 3, 4),
    y = c(3, 4, 6, 8)
  ))
})



test_that("getProblem is available in reduceResults", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) 2)
  addExperiments(reg, "p1", "a1", repls = 2)
  submitJobs(reg)
  waitForJobs(reg)
  z = reduceResults(reg, fun = function(aggr, job, res) {
    pid = job$prob.id
    p = getProblem(reg, pid)
    static = p$static
    list(pid, static, res)
  })
  expect_equal(z, list("p1", 1, 2))
})



test_that("reduceResultsExperiments works with default fun", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) list(y = 1))
  addExperiments(reg, "p1", "a1", repls = 2)
  submitJobs(reg)
  waitForJobs(reg)
  z = reduceResultsExperiments(reg)
  expect_equal(z, data.frame(
    id = 1:2,
    prob = "p1",
    algo = "a1",
    repl = 1:2,
    y = 1), check.attributes = FALSE)

  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) 1)
  addExperiments(reg, "p1", "a1", repls = 2)
  submitJobs(reg)
  waitForJobs(reg)
  z = reduceResultsExperiments(reg)
  expect_equal(z, data.frame(
    id = 1:2,
    prob = "p1",
    algo = "a1",
    repl = 1:2,
    X1 = 1), check.attributes = FALSE)

  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) c(foo = 1, bar = 2))
  addExperiments(reg, "p1", "a1", repls = 2)
  submitJobs(reg)
  waitForJobs(reg)
  z = reduceResultsExperiments(reg, ids = 2)
  data2 = setRowNames(data.frame(
    id = 2L,
    prob = "p1",
    algo = "a1",
    repl = 2L,
    foo = 1,
    bar = 2
  ), 2L)
  attr(data2, "prob.pars.names") = character(0)
  attr(data2, "algo.pars.names") = character(0)
  class(data2) = c("ReducedResultsExperiments", "data.frame")
  expect_equal(z, data2)
})


test_that("reduceResultsExperiments works with imputation", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) list(y = 1))
  ids = addExperiments(reg, "p1", "a1", repls = 3)
  submitJobs(reg, 1:2)
  waitForJobs(reg)
  z = reduceResultsExperiments(reg, ids, impute.val = list())
  expect_equal(z, data.frame(
    id = 1:3,
    prob = "p1",
    algo = "a1",
    repl = 1:3,
    y = c(1, 1, NA)), check.attributes = FALSE)
  z = reduceResultsExperiments(reg, ids, impute.val = list(missing = TRUE))
  expect_equal(z, data.frame(
    id = 1:3,
    prob = "p1",
    algo = "a1",
    repl = 1:3,
    y = c(1, 1, NA),
    missing = c(NA, NA, TRUE)), check.attributes = FALSE)

  # we had a bad bug, where job ids in result df where not correct when one used imputation
  # and the imputed jobs occured before the df row in question
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) list(y = 1))
  ids = addExperiments(reg, "p1", "a1", repls = 3)
  submitJobs(reg, 2:3)
  waitForJobs(reg)
  z = reduceResultsExperiments(reg, ids, impute.val = list())
  expect_equal(z, data.frame(
    id = 1:3,
    prob = "p1",
    algo = "a1",
    repl = 1:3,
    y = c(NA, 1, 1)), check.attributes = FALSE)
  z = reduceResultsExperiments(reg, ids, impute.val = list(missing = TRUE))
  expect_equal(z, data.frame(
    id = 1:3,
    prob = "p1",
    algo = "a1",
    repl = 1:3,
    y = c(NA, 1, 1),
    missing = c(TRUE, NA, NA)), check.attributes = FALSE)
})

test_that("reduceResultsExperiments applies on missing results", {
  reg = makeTestRegistry()
  addProblem(reg, "p1", static = 1)
  addAlgorithm(reg, id = "a1", fun = function(static, dynamic) list(y = 1))
  ids = addExperiments(reg, "p1", "a1", repls = 3)
  fun = function(job, res) list(jid = job$id)
  expect_error(reduceResultsExperiments(reg, ids = 1:3, fun = fun, apply.on.missing = FALSE, "no results"))
  res = reduceResultsExperiments(reg, ids = 1:3, fun = fun, apply.on.missing = TRUE)
  expect_is(res, "data.frame")
  expect_identical(res$jid, 1:3)
})
