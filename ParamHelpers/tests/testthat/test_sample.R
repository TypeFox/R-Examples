context("sample")

test_that("numeric sampling", {
  p = makeNumericParam(id = "x", lower = 10, upper = 20)
  r = sampleValues(p, 13)
  expect_true(is.list(r))
  r = unlist(r)
  expect_true(is.numeric(r))
  expect_equal(length(r), 13)
  expect_true(all(!is.na(r)))
  expect_true(all(r >= p$lower & r <= p$upper))
})

test_that("numeric vector sampling", {
  p = makeNumericVectorParam(id = "x", len = 2, lower = 10, upper = 20)
  r = sampleValues(p, 1000)
  r = do.call(rbind, r)
  expect_true(all(r >= p$lower & r <= p$upper))
})


test_that("integer sampling", {
  p = makeIntegerParam(id = "x", lower = 10, upper = 20)
  r = sampleValues(p, 13)
  expect_true(is.list(r))
  r = unlist(r)
  expect_true(is.integer(r))
  expect_equal(length(r), 13)
  expect_true(all(!is.na(r)))
  expect_true(all(r >= p$lower & r <= p$upper))
})

test_that("integer vector sampling", {
  p = makeIntegerVectorParam(id = "x", len = 2, lower = 1, upper = 3)
  r = sampleValues(p, 1000)
  r = do.call(rbind, r)
  expect_true(all(r >= p$lower & r <= p$upper))
  # this is stochastic, we dont want that on CRAN as it can fail
  if (interactive()) {
    r = as.numeric(table(r))
    expect_true(all(r > 600 & r < 730))
  }
})

test_that("logical sampling", {
  p = makeLogicalParam(id = "x")
  r = sampleValues(p, 13)
  expect_true(is.list(r))
  r = unlist(r)
  expect_true(is.logical(r))
  expect_equal(length(r), 13)
  expect_true(all(!is.na(r)))
})

test_that("logical vector sampling", {
  p = makeLogicalVectorParam(id = "x", len = 2)
  r = sampleValues(p, 1000)
  expect_true(is.list(r))
  expect_true(all(sapply(r, function(x) is.logical(x) && length(x) == 2)))
  expect_true(setequal(unique(unlist(r)), c(TRUE, FALSE)))
  r = do.call(rbind, r)
  if (interactive()) {
    r1 = as.numeric(table(r[,1]))
    expect_true(all(r1 > 300))
    r2 = as.numeric(table(r[,2]))
    expect_true(all(r2 > 300))
  }
})

test_that("discrete sampling", {
  p = makeDiscreteParam(id = "x", values = c("a", "b", "c"))
  r = sampleValues(p, 13)
  expect_true(is.list(r))
  r = unlist(r)
  expect_true(is.character(r))
  expect_equal(length(r), 13)
  expect_true(all(!is.na(r)))
  expect_true(all(r %in% p$values))

  p = makeDiscreteParam(id = "x", values = c(xx = "a", yy = "b"))
  r = sampleValues(p, 13, discrete.names = TRUE)
  expect_true(is.list(r))
  r = unlist(r)
  expect_true(is.character(r))
  expect_equal(length(r), 13)
  expect_true(all(!is.na(r)))
  expect_true(all(r %in% c("xx", "yy")))

  p = makeDiscreteParam(id = "x", values = list(a = list(), b = 1:3))
  r = sampleValues(p, 10)
  expect_true(is.list(r))
  expect_equal(length(r), 10)
  expect_true(all(r %in% p$values))

  p = makeDiscreteVectorParam(id = "x", len = 2, values = list(a = list(), b = 1:3))
  r = sampleValues(p, 10)
  expect_true(is.list(r))
  expect_equal(length(r), 10)
  ok = function(x) is.list(x) && length(x) == 2 &&
    (length(x[[1]]) == 0L || x[[1]] %in% 1:3) &&  (length(x[[2]]) == 0L || x[[2]] %in% 1:3)
  expect_true(all(sapply(r, ok)))
})

test_that("bounds checked", {
  expect_error(sampleValue(makeNumericParam("u", lower = 2)), "Cannot sample")
  expect_error(sampleValue(makeNumericParam("u", upper = 2)), "Cannot sample")
  expect_error(sampleValue(makeIntegerVectorParam("u", len = 2, lower = 2)), "Cannot sample")
  expect_error(sampleValue(makeIntegerVectorParam("u", len = 2, upper = 2)), "Cannot sample")
})

test_that("sampleValues with paramset works", {
  ps = makeParamSet(
    makeIntegerParam("x", lower = 20, upper = 100),
    makeDiscreteParam("y", values = c(xx = "a", yy = "b"))
  )
  r = sampleValues(ps, 5, discrete.names = TRUE)
  expect_true(is.list(r))
  expect_equal(length(r), 5)
  expect_true(all(sapply(r, function(x) is.numeric(x[[1]]) && x[[1]] >= 1 &&  x[[1]] <= 100)))
  expect_true(all(sapply(r, function(x) x[[2]] %in% c("xx", "yy"))))
})


test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower = 1, upper = 2, requires = quote(x == "a")),
    makeIntegerVectorParam("z", len = 2, lower = 1, upper = 20, requires = quote(x == "b")),
    makeDiscreteParam("w", values = 1:2, requires = quote(x == "a"))
  )
  vals = sampleValues(ps, n = 100)
  oks = sapply(vals, function(v) isFeasible(ps, v))
  expect_true(all(oks))
})

test_that("trafo works", {
  ps = makeParamSet(
    makeNumericParam("x", lower =-2, upper = 2, trafo = function(x) 2^x)
  )
  vals = sampleValues(ps, n = 100, trafo = TRUE)
  vals = extractSubList(vals, "x")
  expect_true(all(vals >= 2^(-2) & vals <= 2^2))
})



