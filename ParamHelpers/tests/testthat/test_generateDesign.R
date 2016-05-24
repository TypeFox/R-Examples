context("generateDesign")

test_that("simple num design", {
  requirePackages("_lhs")
  ps1 = makeParamSet(
    makeNumericParam("x1", lower = -2, upper = 1)
  )
  des = generateDesign(13, ps1, lhs::randomLHS)
  expect_equal(nrow(des), 13)
  expect_equal(ncol(des), 1)
  expect_true(is.numeric(des[,1]))
  expect_true(des[,1] >= -2 && des[,1] <= 1)
  des = generateDesign(13, ps1, lhs::maximinLHS)
  expect_equal(nrow(des), 13)
  expect_equal(ncol(des), 1)
  expect_true(des[,1] >= -2 && des[,1] <= 1)
})

test_that("simple num/int design", {
  requirePackages("_lhs")
  ps2 = makeParamSet(
    makeNumericParam("x1", lower = -2, upper = 1),
    makeIntegerParam("x2", lower = 10, upper = 20)
  )
  des = generateDesign(13, ps2, lhs::randomLHS)
  expect_equal(nrow(des), 13)
  expect_equal(ncol(des), 2)
  expect_true(is.numeric(des[,1]))
  expect_true(des[,1] >= -2 && des[,1] <= 1)
  expect_true(is.integer(des[,2]))
  expect_true(des[,2] >= 10 && des[,2] <= 20)
  des = generateDesign(13, ps2, lhs::maximinLHS)
  expect_equal(nrow(des), 13)
  expect_equal(ncol(des), 2)
  expect_true(is.numeric(des[,1]))
  expect_true(des[,1] >= -2 && des[,1] <= 1)
  expect_true(is.integer(des[,2]))
  expect_true(des[,2] >= 10 && des[,2] <= 20)
})

test_that("simple num/int/discrete/log design", {
  ps3 = makeParamSet(
    makeNumericParam("x1", lower = -2, upper = 1),
    makeIntegerParam("x2", lower = 10, upper = 20),
    makeDiscreteParam("x3", values = c("a", "b", "c")),
    makeLogicalParam("x4")
  )
  des = generateDesign(500, ps3)
  expect_equal(nrow(des), 500)
  expect_equal(ncol(des), 4)
  expect_true(is.numeric(des[,1]))
  expect_true(des[,1] >= -2 && des[,1] <= 1)
  expect_true(is.integer(des[,2]))
  expect_true(des[,2] >= 10 && des[,2] <= 20)
  expect_true(is.factor(des[,3]))
  expect_true(all(des[,3] %in% names(ps3$pars[[3]]$values)))
  expect_equal(levels(des[,3]), names(ps3$pars[[3]]$values))
  expect_true(is.logical(des[,4]))
  tab = as.numeric(table(des[,3]))
  expect_true(all(tab > 140 & tab < 180))
  tab = as.numeric(table(des[,4]))
  expect_true(all(tab > 200 & tab < 300))
})

test_that("num/int/disc vec design", {
  ps4 = makeParamSet(
    makeNumericVectorParam("x", len = 2, lower = -2, upper = 1),
    makeIntegerVectorParam("y", len = 3, lower = 10L, upper = 20L),
    makeDiscreteVectorParam("z", len = 2, values = list(a = "a", b = list())),
    makeLogicalVectorParam("a", len = 2)
  )
  des = generateDesign(13, ps4)
  expect_equal(nrow(des), 13)
  expect_equal(ncol(des), 9)
  expect_equal(colnames(des), c("x1", "x2", "y1", "y2", "y3", "z1", "z2", "a1", "a2"))
  expect_true(is.numeric(des[,1]))
  expect_true(is.numeric(des[,2]))
  expect_true(is.integer(des[,3]))
  expect_true(is.integer(des[,4]))
  expect_true(is.integer(des[,5]))
  expect_true(is.factor(des[,6]))
  expect_true(is.factor(des[,7]))
  expect_true(is.logical(des[,8]))
  expect_true(is.logical(des[,9]))
  expect_true(des[,1] >= -2 && des[,1] <= 1)
  expect_true(des[,2] >= -2 && des[,2] <= 1)
  expect_true(des[,3] >= 10 && des[,3] <= 20)
  expect_true(des[,4] >= 10 && des[,4] <= 20)
  expect_true(des[,5] >= 10 && des[,5] <= 20)
  expect_true(all(des[,6] %in% c("a", "b")))
  expect_true(all(des[,7] %in% c("a", "b")))
  expect_equal(levels(des[,6]), c("a", "b"))
  expect_equal(levels(des[,7]), c("a", "b"))
})

test_that("num/int vec design with trafo", {
  ps5 = makeParamSet(
    makeNumericVectorParam("x", len = 2, lower = -2, upper = 1, trafo =  function(x) 2^x),
    makeIntegerVectorParam("y", len = 3, lower = 10L, upper = 20L, trafo = function(x) 3L*x)
  )
  des = generateDesign(100, ps5, trafo = TRUE)
  expect_equal(nrow(des), 100)
  expect_equal(ncol(des), 5)
  expect_equal(colnames(des), c("x1", "x2", "y1", "y2", "y3"))
  expect_true(is.numeric(des[,1]))
  expect_true(is.numeric(des[,2]))
  expect_true(is.integer(des[,3]))
  expect_true(is.integer(des[,4]))
  expect_true(is.integer(des[,5]))
  expect_true(des[,1] >= 1/4 && des[,1] <= 2)
  expect_true(des[,2] >= 1/4 && des[,2] <= 2)
  expect_true(des[,3] >= 30 && des[,3] <= 60)
  expect_true(des[,4] >= 30 && des[,4] <= 60)
  expect_true(des[,5] >= 30 && des[,5] <= 60)

  ps6 = makeParamSet(
    makeNumericVectorParam("x", len = 2, lower = 0, upper = 1, trafo = function(x) x / sum(x)),
    makeIntegerVectorParam("y", len = 2, lower = 10L, upper = 20L, trafo = function(x) 1:2),
    makeDiscreteVectorParam("z", len = 2, values = c("a", "b"))
  )
  des = generateDesign(5, ps6, trafo = TRUE)
  expect_equal(nrow(des), 5)
  expect_equal(ncol(des), 6)
  expect_equal(colnames(des), c("x1", "x2", "y1", "y2", "z1", "z2"))
  expect_true(is.numeric(des[,1]))
  expect_true(is.numeric(des[,2]))
  expect_true(is.integer(des[,3]))
  expect_true(is.integer(des[,4]))
  expect_equal(des[,1] + des[,2], rep(1, 5))
  expect_equal(des[,3], rep(1, 5))
  expect_equal(des[,4], rep(2, 5))

  # check that trafo can be switched off
  ps = makeParamSet(
    makeNumericParam("x", lower = -10, upper = -1, trafo = function(x) -x)
  )
  des = generateDesign(10, ps, trafo = FALSE)
  expect_true(is.numeric(des$x) && all(des$x <= -1))
})

test_that("requires works", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericParam("y", lower = 1, upper = 2, requires = quote(x == "a")),
    makeDiscreteParam("z", values = 1:2, requires = quote(x == "b"))
  )
  des = generateDesign(50, par.set = ps)
  vals = dfRowsToList(des, ps)
  oks = sapply(vals, isFeasible, par = ps)
  expect_true(all(oks))
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b")),
    makeNumericVectorParam("y", len = 2, lower = 1, upper = 2, requires = quote(x == "a"))
  )
  des = generateDesign(50, par.set = ps)
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
  des = generateDesign(50, par.set = ps7)
  expect_true(all(is.na(des[des$disc == "a", 5:8])))
  expect_true(all(is.na(des[des$disc == "b", 2:4])))
  expect_true(all(is.na(des[des$disc == "c", 2:8])))
  expect_true(all(is.na(des[des$disc == "b" & des$discB == "NR",7])))
  expect_true(all(is.na(des[des$disc == "b" & des$discB == "R",8])))

  vals = dfRowsToList(des, ps7)
  oks = sapply(vals, isFeasible, par = ps7)
  expect_true(all(oks))
})

test_that("requires chains work", {
  ps8 = makeParamSet(
    makeLogicalLearnerParam("a", default = FALSE),
    makeLogicalLearnerParam("b", default = FALSE, requires = quote(a == TRUE)),
    # FIXME: shouldn't need to make the chain explicit
    makeLogicalLearnerParam("c", default = FALSE, requires = quote(b == TRUE && a == TRUE))
  )
  des = generateDesign(10, par.set = ps8)
  vals = dfRowsToList(des, ps8)
  oks = sapply(vals, isFeasible, par = ps8)
  expect_true(all(oks))
})


test_that("we dont drop levels in factors", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = letters[5:1])
  )
  des = generateDesign(1, ps)
  expect_true(is.factor(des$x) && levels(des$x) == letters[5:1])
})

test_that("list of further arguments passed to lhs function produces no error", {
  ps = makeParamSet(
    makeDiscreteParam("x", values = c("a", "b", "c")),
    makeNumericParam("realA", lower = 0, upper = 100)
  )
  generateDesign(10L, ps, fun.args = list(preserveDraw = TRUE))
})

test_that("removing duplicates", {
  # this param set has 6 possilbe combinations of param values ...
  ps = makeParamSet(
    makeIntegerParam("int", lower = 0, upper = 1),
    makeDiscreteParam("disc", values = letters[1:3])
  )
  # ... and therefore creating a design with n > 6 duplicate free rows should fail
  expect_warning(res <- generateDesign(7L, ps))
  expect_true(nrow(res) <= 6L)
})
