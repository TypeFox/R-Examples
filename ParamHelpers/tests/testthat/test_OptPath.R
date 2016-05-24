context("OptPath")

test_that("OptPath", {
  ps = makeParamSet(
    makeNumericParam("x"),
    makeDiscreteParam("y", values = c("a", "b"))
  )
  op = makeOptPathDF(par.set = ps, y.names = c("z1", "z2"), minimize = c(TRUE, FALSE))
  addOptPathEl(op, x = list(x = 1, y = "a"), y = c(z1 = 1, z2 = 4))
  addOptPathEl(op, x = list(x = 2, y = "a"), y = c(z1 = 3, z2 = 2))
  expect_equal(op$env$dob, 1:2)
  setOptPathElEOL(op, 2, 8)
  expect_equal(op$env$eol[2], 8)

  # test getters
  expect_equal(getOptPathX(op), data.frame(x = 1:2, y = "a"))
  expect_equal(getOptPathX(op, dob = 2), data.frame(x = 2, y = "a"))

  expect_equal(getOptPathY(op, "z1"), c(1, 3))
  expect_equal(getOptPathY(op, "z2"), c(4, 2))
  expect_equal(getOptPathY(op), matrix(c(1, 3, 4, 2), nrow = 2L, dimnames = list(1:2, c("z1", "z2"))))
  expect_equal(getOptPathY(op, "z2", drop = FALSE), matrix(c(4, 2), nrow = 2L, dimnames = list(1:2, c("z2"))))
  expect_equal(getOptPathY(op, "z2", drop = TRUE), c(4, 2))
  expect_equal(getOptPathY(op, "z2"), c(4, 2))
  expect_equal(getOptPathY(op, eol = 8), matrix(c(3, 2), nrow = 1L, dimnames = list(2, c("z1", "z2"))))
  expect_equal(getOptPathDOB(op), c(1, 2))
  expect_equal(getOptPathEOL(op), c(NA, 8))

  x = as.data.frame(op)
  expect_true(is.data.frame(x))
  expect_equal(nrow(x), 2)
  expect_equal(ncol(x), 6)

  expect_output(print(op), "Optimization path")

  expect_equal(getOptPathEl(op, 1)$x, list(x = 1, y = "a"))

  gbe = function(op, y.name, dob) {
    i = getOptPathBestIndex(op, y.name, dob)
    getOptPathEl(op, i)
  }

  expect_equal(gbe(op, y.name = "z1", dob = 1:2), getOptPathEl(op, 1))
  expect_equal(gbe(op, y.name = "z2", dob = 1:2), getOptPathEl(op, 1))
  expect_equal(gbe(op, y.name = "z1", dob = 1), getOptPathEl(op, 1))
  expect_equal(gbe(op, y.name = "z2", dob = 1), getOptPathEl(op, 1))
  expect_equal(gbe(op, y.name = "z1", dob = 2), getOptPathEl(op, 2))
  expect_equal(gbe(op, y.name = "z2", dob = 2), getOptPathEl(op, 2))

  setOptPathElDOB(op, 1, 1)
  setOptPathElDOB(op, 2, 3)
  expect_equal(as.data.frame(op)$dob, c(1, 3))
  setOptPathElDOB(op, 1:2, 7)
  expect_equal(as.data.frame(op)$dob, c(7, 7))
  setOptPathElDOB(op, 1:2, 4:5)
  expect_equal(as.data.frame(op)$dob, c(4, 5))
  setOptPathElEOL(op, 1, 1)
  setOptPathElEOL(op, 2, 3)
  expect_equal(as.data.frame(op)$eol, c(1, 3))
  setOptPathElEOL(op, 1:2, 7)
  expect_equal(as.data.frame(op)$eol, c(7, 7))
  setOptPathElEOL(op, 1:2, 4:5)
  expect_equal(as.data.frame(op)$eol, c(4, 5))

  ps = makeParamSet(
    makeNumericVectorParam("x", len = 2),
    makeIntegerParam("y")
  )
  op = makeOptPathDF(par.set = ps, y.names = "z", minimize = TRUE)
  addOptPathEl(op, x = list(c(1,1), 7L), y = 1)
  addOptPathEl(op, x = list(c(2,2), 8L), y = 3)
  df = as.data.frame(op)
  expect_equal(dim(df), c(2, 3 + 1 + 2))
  expect_equal(colnames(df), c("x1", "x2", "y", "z", "dob", "eol"))
  e = getOptPathEl(op, 1)
  expect_equal(e$x, list(x = c(1,1), y = 7L))
  # really make sure that names are there
  expect_equal(names(e$y), "z")
  # check error msg for wrong y
  expect_error(addOptPathEl(op, x = list(c(1,1), 7L), y = c(1, 1)), "Must have length")
})

test_that("OptPath with vector and discrete params works", {
  ps = makeParamSet(
    makeIntegerVectorParam("x0", len = 2),
    makeDiscreteParam("x1", values = c("a", "b")),
    makeDiscreteParam("x2", values = 1:2),
    makeDiscreteParam("x3", values = c(1.2, 5)),
    makeDiscreteParam("x4", values = list(foo = identity, bar = list())),
    makeLogicalParam("x5"),
    makeDiscreteVectorParam("x6", len = 2, values = list(a = identity, b = 1)),
    makeLogicalVectorParam("x7", len = 3),
    makeCharacterParam("x8")
  )
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE)
  x = list(x0 = c(1L, 3L), x1 = "a", x2 = 2L, x3 = 5, x4 = identity, x5 = FALSE,
    x6 = list(b = 1, a = identity), x7 = c(TRUE, FALSE, TRUE), x8 = "PH")
  addOptPathEl(op, x = x, y = 0)
  d = as.data.frame(op)
  expect_true(nrow(d) == 1 && ncol(d) == 3 + 8 + 1 + 2 + 2)
  expect_true(is.integer(d$x01))
  expect_true(is.integer(d$x02))
  expect_true(is.factor(d$x1))
  expect_true(is.factor(d$x2))
  expect_true(is.factor(d$x3))
  expect_true(is.factor(d$x4))
  expect_true(is.logical(d$x5))
  expect_true(is.factor(d$x61))
  expect_true(is.factor(d$x62))
  expect_true(is.logical(d$x71))
  expect_true(is.logical(d$x72))
  expect_true(is.logical(d$x73))
  expect_true(is.factor(d$x8)) # strings internally saved as factors?
  expect_true(
       d[1,1] == 1L && d[1,2] == 3L
    && d[1,3] == "a" && d[1,4] == "2" && d[1,5] == "5" && d[1,6] == "foo"
    && d[1,7] == FALSE && d[1,8] == "b" && d[1,9] == "a"
    && d[1,10] == TRUE && d[1,11] == FALSE && d[1,12] == TRUE
    && d[1,13] == "PH")
  d = as.data.frame(op, discretes.as.factor = TRUE)
  expect_true(nrow(d) == 1 && ncol(d) == 3 + 8 + 1 + 2 + 2)
  expect_true(is.integer(d$x01))
  expect_true(is.integer(d$x02))
  expect_true(is.factor(d$x1))
  expect_true(is.factor(d$x2))
  expect_true(is.factor(d$x3))
  expect_true(is.factor(d$x4))
  expect_true(is.logical(d$x5))
  expect_true(is.factor(d$x61))
  expect_true(is.factor(d$x62))
  expect_true(is.logical(d$x71))
  expect_true(is.logical(d$x72))
  expect_true(is.logical(d$x73))
  expect_true(is.factor(d$x8)) # strings internally saved as factors?
  expect_true(
    d[1,1] == 1L && d[1,2] == 3L
    && d[1,3] == "a" && d[1,4] == "2" && d[1,5] == "5" && d[1,6] == "foo"
    && d[1,7] == FALSE && d[1,8] == "b" && d[1,9] == "a"
    && d[1,10] == TRUE && d[1,11] == FALSE && d[1,12] == TRUE
    && d[1,13] == "PH")
  e = getOptPathEl(op, 1L)
  expect_equal(e$x, x)
})

test_that("getOptPathBestIndex tie-handling", {
  ps = makeParamSet(
    makeNumericParam("x")
  )
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 1)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 1)
  expect_equal(getOptPathBestIndex(op), 4L)
  expect_equal(getOptPathBestIndex(op, ties = "first"), 1L)
  expect_equal(getOptPathBestIndex(op, ties = "last"), 4L)
  expect_true(getOptPathBestIndex(op, ties = "random") %in% c(1L, 2L, 4L))
  expect_equal(getOptPathBestIndex(op, ties = "all"), c(1L, 2L, 4L))

  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = FALSE)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 1)
  addOptPathEl(op, x = list(x = 0), y = 0)
  addOptPathEl(op, x = list(x = 0), y = 1)
  expect_equal(getOptPathBestIndex(op), 5L)
  expect_equal(getOptPathBestIndex(op, ties = "first"), 3L)
  expect_equal(getOptPathBestIndex(op, ties = "last"), 5L)
  expect_true(getOptPathBestIndex(op, ties = "random") %in% c(3L, 5L))
  expect_equal(getOptPathBestIndex(op, ties = "all"), c(3L, 5L))
})


test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower = 1, upper = 2, requires = quote(x == "a")),
    makeIntegerVectorParam("z", len = 2, lower = 1, upper = 20, requires = quote(x == "b"))
  )
  op = makeOptPathDF(par.set = ps, y.names = "foo", minimize = TRUE)
  el = list(x = "b", y = NA, z = 1:2)
  addOptPathEl(op, x = el, y = 0)
  # check that correct NAs are present
  expect_equal(op$env$path[, 2:4], data.frame(y = NA_real_, z1 = 1L, z2 = 2L))

  op = makeOptPathDF(par.set = ps, y.names = "foo", minimize = TRUE)
  el = list(x = "a", y = 1, z = NA)
  addOptPathEl(op, x = el, y = 0)
  expect_equal(getOptPathEl(op, 1)$x, el)
  # check that correct NAs are present
  expect_equal(op$env$path[, 2:4], data.frame(y = 1, z1 = NA_integer_, z2 = NA_integer_))
  el = list(x = "b", y = NA, z = c(2, 3))
  addOptPathEl(op, x = el, y = 0)
  expect_equal(getOptPathEl(op, 2)$x, el)
  d = as.data.frame(op, discretes.as.factor = TRUE)
  expect_equal(d[,1:4], data.frame(x = c("a", "b"), y = c(1, NA), z1 = c(NA, 2L), z2 = c(NA, 3L)))
})

test_that("pareto front", {
  ps = makeParamSet(makeNumericParam("x"))
  op = makeOptPathDF(par.set = ps, y.names = c("y1", "y2"), minimize = c(TRUE, TRUE))
  addOptPathEl(op, x = list(x = 3), y = c(9, 4))
  f = getOptPathParetoFront(op, index = FALSE)
  expect_equal(f, matrix(c(9, 4), nrow = 1L, dimnames = list(1, c("y1", "y2"))))
  addOptPathEl(op, x = list(x = 4), y = c(4, 9))
  addOptPathEl(op, x = list(x = 1), y = c(5, 3))
  addOptPathEl(op, x = list(x = 2), y = c(2, 4))
  f = getOptPathParetoFront(op, index = TRUE)
  expect_equal(f, 3:4)
  f = getOptPathParetoFront(op, index = FALSE)
  expect_equal(f, matrix(c(5, 2, 3, 4), nrow = 2L, dimnames = list(3:4, c("y1", "y2"))))
  f = getOptPathParetoFront(op, index = TRUE, dob = 3:4)
  expect_equal(f, 3:4)
  f = getOptPathParetoFront(op, index = FALSE, dob = 3:4)
  expect_equal(f, matrix(c(5, 2, 3, 4), nrow = 2L, dimnames = list(3:4, c("y1", "y2"))))
  op$minimize = c(TRUE, FALSE)
  f = getOptPathParetoFront(op, index = TRUE)
  expect_equal(f, c(2, 4))
})

test_that("error message and exec time works", {
  ps = makeParamSet(makeNumericParam("x"))
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE, include.error.message = TRUE)
  addOptPathEl(op, x = list(x = 1), y = 5)
  addOptPathEl(op, x = list(x = 2), y = 3)
  addOptPathEl(op, x = list(x = 3), y = 9)
  addOptPathEl(op, x = list(x = 4), y = 3, error.message = "bla")
  expect_equal(ncol(as.data.frame(op)), 5)
  errors = which(!is.na(getOptPathErrorMessages(op)))
  expect_equal("bla", getOptPathEl(op, errors)$error.message)

  ps = makeParamSet(makeNumericParam("x"))
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE, include.error.message = TRUE,
    include.exec.time = TRUE)
  addOptPathEl(op, x = list(x = 1), y = 5, exec.time = 2)
  addOptPathEl(op, x = list(x = 2), y = 3, exec.time = 2)
  addOptPathEl(op, x = list(x = 3), y = 9, exec.time = 3)
  addOptPathEl(op, x = list(x = 4), y = 3, exec.time = NA_real_, error.message = "bla")
  expect_equal(ncol(as.data.frame(op)), 6)
  errors = which(!is.na(getOptPathErrorMessages(op)))
  expect_equal("bla", getOptPathEl(op, errors)$error.message)
  expect_equal(getOptPathExecTimes(op), c(2, 2, 3, NA))
})

test_that("logging extra works", {
  ps = makeParamSet(makeNumericParam("v"))
  op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE, include.extra = TRUE)
  addOptPathEl(op, x = list(v = 1), y = 5, extra = list(ee = 7))
  df = setRowNames(as.data.frame(op), NULL)
  expect_equal(df, data.frame(v = 1, y = 5, dob = 1L, eol = NA_integer_, ee = 7))
  expect_equal(getOptPathEl(op, 1L), list(x = list(v = 1), y = c(y = 5), dob = 1L, eol = NA_integer_,
    extra = list(ee = 7)))
})

test_that("as.data.frame flags and getCols works", {
  ps = makeParamSet(
    makeNumericParam("x"),
    makeDiscreteParam("y", values = c("a", "b"))
  )
  op = makeOptPathDF(par.set = ps, y.names = c("z1", "z2"), minimize = c(TRUE, TRUE), include.extra = TRUE)
  addOptPathEl(op, x = list(x = 1, y = "a"), y = c(z1 = 1, z2 = 4), extra = list(ee = 7))
  addOptPathEl(op, x = list(x = 2, y = "a"), y = c(z1 = 3, z2 = 2), extra = list(ee = 8))

  expect_error(as.data.frame(op, include.x = FALSE, include.y = FALSE, include.rest = FALSE), "include something")
  df1 = as.data.frame(op, include.rest = FALSE, discretes.as.factor = TRUE)
  df2 = as.data.frame(op, include.rest = FALSE, include.x = FALSE, discretes.as.factor = TRUE)
  df3 = as.data.frame(op, include.y = FALSE, include.x = FALSE, discretes.as.factor = TRUE)
  df4 = as.data.frame(op, include.y = TRUE, include.x = TRUE, include.rest = TRUE, dob = 2L)
  expect_equal(nrow(df1), 2)
  expect_equal(nrow(df2), 2)
  expect_equal(nrow(df3), 2)
  expect_equal(nrow(df4), 1)
  expect_equal(ncol(df2), 2)
  expect_equal(ncol(df3), 3)
  expect_equal(ncol(df4), 7)
  expect_equal(sapply(df1, class), c(x = "numeric", y = "factor", z1 = "numeric", z2 = "numeric"))
  expect_equal(sapply(df2, class), c(z1 = "numeric", z2 = "numeric"))
  expect_equal(sapply(df3, class), c(dob = "integer", eol = "integer", ee = "numeric"))
  expect_equal(sapply(df4, class), c(x = "numeric", y = "factor", z1 = "numeric", z2 = "numeric",
    dob = "integer", eol = "integer", ee = "numeric"))
  expect_error(getOptPathCol(op,"bla"), "not present")
  expect_equal(getOptPathCol(op, "dob"), c(1, 2))
  expect_equal(getOptPathCol(op, "eol"), c(NA_integer_, NA_integer_))
  expect_equal(getOptPathCol(op, "error.message"), NULL)
  expect_equal(getOptPathCol(op, "exec.time"), NULL)
  expect_equal(getOptPathCol(op, "x"), 1:2)
  expect_equal(getOptPathCol(op, "y"), c("a", "a"))
  expect_equal(getOptPathCol(op, "z1"), c(1, 3))
  expect_equal(getOptPathCol(op, "ee"), c(7, 8))

  d = getOptPathCols(op, c("x", "y"))
  expect_equal(dim(d), c(2L, 2L))
  expect_equal(colnames(d), c("x", "y"))
  expect_true(is.factor(d$y))
  d = getOptPathCols(op, "y")
  expect_equal(dim(d), c(2L, 1L))
  expect_true(is.factor(d$y))

})


test_that("opt.path printing works", {
  ps = makeParamSet(
    makeNumericParam("x"),
    makeDiscreteParam("y", values = c("a", "b"))
  )
  op = makeOptPathDF(par.set = ps, y.names = c("z1", "z2"), minimize = c(TRUE, FALSE))
  expect_output(print(op), "Optimization path")
  addOptPathEl(op, x = list(x = 1, y = "a"), y = c(z1 = 1, z2 = 4))
  expect_output(print(op), "Optimization path")

  op = makeOptPathDF(par.set = ps, y.names = "z1", minimize = TRUE, include.exec.time = TRUE)
  expect_output(print(op), "Exec times: TRUE. Range: 0 - 0. 0 NAs")
  addOptPathEl(op, x = list(x = 1, y = "a"), y = c(z1 = 1))
  expect_output(print(op), "Exec times: TRUE. Range: 0 - 0. 1 NAs")
  addOptPathEl(op, x = list(x = 1, y = "a"), y = c(z1 = 1), exec.time = 3)
  expect_output(print(op), "Exec times: TRUE. Range: 3 - 3. 1 NAs.")
})
