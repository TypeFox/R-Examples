context("getJobInfo")

test_that("getJobInfo", {
  r = makeTestRegistry()
  p1 = addProblem(r, "one", 1)
  p2 = addProblem(r, "two", 2)
  a1 = addAlgorithm(r, "A", fun=function(static, dynamic) 1)
  a2 = addAlgorithm(r, "B", fun=function(static, dynamic) 1)
  addExperiments(r, list(makeDesign(p1), makeDesign(p2)), list(makeDesign(a1), makeDesign(a2)), repls=2)

  mycheck = function(tab) {
    expect_true(is.data.frame(tab))
    expect_equal(tab$id, 1:8)
    expect_true(nrow(tab) == 8)
    expect_true(is(tab$time.submitted, "POSIXt"))
    expect_true(is(tab$time.started, "POSIXt"))
    expect_true(is(tab$time.done, "POSIXt"))
    expect_true(is.numeric(tab$time.queued))
    expect_true(is.numeric(tab$time.running))
    expect_true(all(is.na(tab$error.msg)))
    expect_true(is.integer(tab$r.pid))
    expect_true(is.integer(tab$seed))
  }
  tab = getJobInfo(r)
  mycheck(tab)
  submitJobs(r)
  waitForJobs(r)
  tab = getJobInfo(r)
  mycheck(tab)

  tab = getJobInfo(r, ids = integer(0))
  expect_true(is.data.frame(tab))
  expect_true(nrow(tab) == 0L)

  tab = getJobInfo(r, ids = 1, pars=TRUE)
  expect_equal(tab$prob, "one")
  expect_equal(tab$algo, "A")

  tab = getJobInfo(r, select = "time.running")
  expect_true(ncol(tab) == 2) # job.id always selected
  tab = getJobInfo(r, select = c("id", "time.running"))
  expect_true(ncol(tab) == 2)

  expect_equal(getJobInfo(r)$repl, rep(1:2, 4))
})
