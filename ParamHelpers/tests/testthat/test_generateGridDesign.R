context("generateGridDesign")

test_that("generateGridDesign", {
  ps = makeParamSet(
    makeNumericParam("x", lower = 1, upper = 5),
    makeIntegerParam("y", lower = 2, upper= 6)
  )
  d = generateGridDesign(ps, resolution = 3L)
  e = expand.grid(x = c(1, 3, 5), y = c(2L, 4L, 6L), KEEP.OUT.ATTRS = FALSE)
  attr(e, "trafo") = FALSE
  expect_equal(d, e)

  ps = makeParamSet(
    makeNumericParam("u", lower = 1, upper = 5),
    makeIntegerParam("v", lower = 2, upper= 6),
    makeLogicalParam("w"),
    makeDiscreteParam("x", values = c("a", "b"))
  )
  d = generateGridDesign(ps, resolution = 3L)
  e = expand.grid(u = c(1, 3, 5), v = c(2L, 4L, 6L), w = c(TRUE, FALSE), x = c("a", "b"),
    KEEP.OUT.ATTRS = FALSE)
  attr(e, "trafo") = FALSE
  expect_equal(d, e)

  # vectors
  ps = makeParamSet(
    makeNumericVectorParam("x", len = 2L, lower = 1, upper = 2),
    makeIntegerVectorParam("y", len = 2L, lower = 3, upper = 4),
    makeLogicalVectorParam("z", len = 2L)
  )
  d = generateGridDesign(ps, resolution = 2L)
  e = expand.grid(
    x1 = c(1, 2), x2 = c(1, 2),
    y1 = c(3L, 4L), y2 = c(3L, 4L),
    z1 = c(TRUE, FALSE), z2 = c(TRUE, FALSE),
    KEEP.OUT.ATTRS = FALSE
  )
  attr(e, "trafo") = FALSE
  expect_equal(d, e)

  # trafo
  ps = makeParamSet(
    makeNumericParam("x", lower = 0, upper = 1),
    makeNumericParam("y", lower = 3, upper = 4, trafo = function(x) 2*x)
  )
  d = generateGridDesign(ps, resolution = c(y = 4, x = 2), trafo = TRUE)
  e = expand.grid(
    x = seq(0, 1, length.out = 2),
    y = 2*(seq(3, 4, length.out = 4)),
    KEEP.OUT.ATTRS = FALSE
  )
  attr(e, "trafo") = TRUE
  expect_equal(d, e)
  ps = makeParamSet(
    makeNumericVectorParam("x", len = 2, lower = 1, upper = 2, trafo = function(x) x / sum(x))
  )
  d = generateGridDesign(ps, resolution = 4, trafo = TRUE)
  e = expand.grid(
    x1 = seq(1, 2, length.out = 4),
    x2 = seq(1, 2, length.out = 4),
    KEEP.OUT.ATTRS = FALSE
  )
  f = as.data.frame(t(apply(e, 1, function(x) x / sum(x))))
  f = f[!duplicated(f),]
  attr(f, "trafo") = TRUE
  expect_equal(setRowNames(d, NULL), setRowNames(f, NULL))
})

test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower = 1, upper = 2, requires = quote(x == "a"))
  )
  des = generateGridDesign(par.set = ps, resolution = 3)
  vals = dfRowsToList(des, ps)
  oks = sapply(vals, isFeasible, par = ps)
  expect_true(all(oks))
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericVectorParam("y", len = 2, lower = 1, upper = 2, requires = quote(x == "a"))
  )
  des = generateGridDesign(par.set = ps, resolution = 3)
  vals = dfRowsToList(des, ps)
  oks = sapply(vals, isFeasible, par = ps)
  expect_true(all(oks))
})

test_that("nested requires", {
  ps7 = makeParamSet(
    makeDiscreteParam("disc", values = c("a", "b", "c")),
    makeNumericParam("realA", lower = 0, upper = 100, requires = quote(disc == "a")),
    makeIntegerParam("intA", lower = -100, upper = 100, requires = quote(disc == "a")),
    makeDiscreteParam("discA", values = c("m", "w"), requires = quote(disc == "a")),
    makeNumericParam("realB", lower = -100, upper = 100, requires = quote(disc == "b")),
    makeDiscreteParam("discB", values = c("R", "NR"), requires = quote(disc == "b")),
    makeNumericParam("realBR", lower = 0, upper = 2*pi, requires = quote(identical(discB, "R") && identical(disc, "b"))),
    makeNumericParam("realBNR", lower = 0, upper = 2*pi, requires = quote(identical(discB, "NR") && identical(disc, "b")))
  )
  des = generateGridDesign(par.set = ps7, resolution = 3)
  expect_true(all(is.na(des[des$disc == "a",5:8])))
  expect_true(all(is.na(des[des$disc == "b",2:4])))
  expect_true(all(is.na(des[des$disc == "c",2:8])))
  expect_true(all(is.na(des[des$disc == "b" & des$discB == "NR",7])))
  expect_true(all(is.na(des[des$disc == "b" & des$discB == "R",8])))

  vals = dfRowsToList(des, ps7)
  oks = sapply(vals, isFeasible, par = ps7)
  expect_true(all(oks))
})

test_that("discrete works without resolution", {
  ps8 = makeParamSet(
    makeDiscreteParam("disc", values = c("c", "b", "a")),
    makeDiscreteParam("discA", values = c("m", "w"), requires = quote(disc == "a")),
    makeLogicalParam("logA")
    )
  des = generateGridDesign(par.set = ps8)
  expect_true(nrow(des) == 8)
  expect_true(all(is.na(des[des$disc == "c", "discA"])))
  expect_true(all(is.na(des[des$disc == "b", "discA"])))
  expect_equal(levels(des[,2]), c("m", "w"))
  expect_equal(levels(des[,1]), c("c", "b", "a"))
})
