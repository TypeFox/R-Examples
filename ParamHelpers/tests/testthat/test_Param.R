context("Parameter")

test_that("num param", {
  p = makeNumericParam(id = "x", lower = -1L, upper = 1)
  expect_equal("numeric", p$type)
  expect_equal(-1, p$lower)
  expect_equal(1, p$upper)

  expect_true(isFeasible(p, -1))
  expect_true(isFeasible(p, 1))
  expect_true(isFeasible(p, 0))

  expect_true(!isFeasible(p, 2))
  expect_true(!isFeasible(p, Inf))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, NA))
  expect_true(!isFeasible(p, "bam"))

  p = makeNumericParam(id = "x", lower = 0, upper = Inf)
  expect_true(isFeasible(p, 2))
  expect_true(!isFeasible(p, -2))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, NULL))

  expect_equal(p$values, NULL)

  # defaults
  p = makeNumericParam(id = "x", allow.inf = TRUE, default = Inf)
  expect_error(makeNumericParam(id = "x", allow.inf = FALSE, default = Inf), "feasible")

  ## Error conditions:
  expect_error(makeNumericParam(id = "x", lower = "bam", upper = 1))
  expect_error(makeNumericParam(id = "x", lower = NA, upper = 1))
  expect_error(makeNumericParam(id = "x", lower = NULL, upper = 1))
  expect_error(makeNumericParam(id = "x", lower = 0, upper = "bam"))
  expect_error(makeNumericParam(id = "x", lower = 0, upper = NA))
  expect_error(makeNumericParam(id = "x", lower = 0, upper = NULL))
  expect_error(makeNumericParam(id = "x", lower = 1, upper = -1))
  expect_error(makeNumericParam(id = "x", lower = c(-1, 1), upper = 2))
  expect_error(makeNumericParam(id = "x", lower = -1, upper = c(1, 2)))
})

test_that("num vec param", {
  p = makeNumericVectorParam(id = "x", lower = -1L, upper = 1, len = 2)
  expect_equal("numericvector", p$type)
  expect_equal(c(-1,-1), p$lower)
  expect_equal(c(1,1), p$upper)

  expect_true(isFeasible(p, c(-1,-1)))
  expect_true(isFeasible(p, c(1,1)))
  expect_true(isFeasible(p, c(0,1)))

  expect_true(!isFeasible(p, c(2,0)))
  expect_true(!isFeasible(p, c(NA,0)))
  expect_true(!isFeasible(p, Inf))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, NA))
  expect_true(!isFeasible(p, "bam"))

  p = makeNumericVectorParam(id = "x", lower = 0, upper = Inf, len = 3)
  expect_true(isFeasible(p, c(2,1,1)))
  expect_true(!isFeasible(p, c(-2,1,0)))
  expect_true(!isFeasible(p, c(-Inf, 1)))
  expect_true(!isFeasible(p, NULL))

  ## Error conditions:
  expect_error(makeNumericVectorParam(id = "x", lower = "bam", upper = 1))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = NA, upper = 1))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = NULL, upper = 1))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = 0, upper = "bam"))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = 0, upper = NA))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = 0, upper = NULL))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = 1, upper = -1))
  expect_error(makeNumericVectorParam(id = "x", len = 2, lower = c(1, 1), upper = c(0,0)))
})

test_that("int param", {
  p = makeIntegerParam(id = "x", lower = -1L, upper = 1)
  expect_equal("integer", p$type)
  expect_equal(-1, p$lower)
  expect_equal(1, p$upper)

  expect_true(isFeasible(p, -1))
  expect_true(isFeasible(p, -1L))
  expect_true(isFeasible(p, 1L))
  expect_true(isFeasible(p, 0L))

  expect_true(!isFeasible(p, 0.5))
  expect_true(!isFeasible(p, Inf))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, NA))
  expect_true(!isFeasible(p, "bam"))

  p = makeIntegerParam(id = "x", lower = 0)
  expect_true(isFeasible(p, 2L))
  expect_true(!isFeasible(p, Inf))
  expect_true(!isFeasible(p, -2L))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, NULL))

  ## Error conditions:
  expect_error(makeIntegerParam(id = "x", lower = "bam", upper = 1L))
  expect_error(makeIntegerParam(id = "x", lower = NA, upper = 1L))
  expect_error(makeIntegerParam(id = "x", lower = NULL, upper = 1L))
  expect_error(makeIntegerParam(id = "x", lower = 0L, upper = "bam"))
  expect_error(makeIntegerParam(id = "x", lower = 0L, upper = NA))
  expect_error(makeIntegerParam(id = "x", lower = 0L, upper = NULL))
  expect_error(makeIntegerParam(id = "x", lower = 1L, upper = -1L))
  expect_error(makeIntegerParam(id = "x", lower = c(-1L, 1L), upper = 2L))
  expect_error(makeIntegerParam(id = "x", lower = -1L, upper = c(1L, 2L)))
})

test_that("int vec param", {
  p = makeIntegerVectorParam(id = "x", lower = -10L, upper = 10, len = 2)
  expect_equal("integervector", p$type)
  expect_equal(c(-10,-10), p$lower)
  expect_equal(c(10,10), p$upper)

  expect_true(isFeasible(p, c(-10,-1)))
  expect_true(isFeasible(p, c(1,1)))
  expect_true(isFeasible(p, c(0,10)))

  expect_true(!isFeasible(p, c(20,0)))
  expect_true(!isFeasible(p, Inf))
  expect_true(!isFeasible(p, -Inf))
  expect_true(!isFeasible(p, c(NA, 5)))
  expect_true(!isFeasible(p, "bam"))

  p = makeIntegerVectorParam(id = "x", lower = 0, len = 3)
  expect_true(isFeasible(p, c(2,1,1)))
  expect_true(isFeasible(p, c(10^7, 10^7, 10^7)))
  expect_true(!isFeasible(p, c(-2,1,0)))
  expect_true(!isFeasible(p, c(-Inf, 1)))
  expect_true(!isFeasible(p, NULL))

  ## Error conditions:
  expect_error(makeIntegerVectorParam(id = "x", lower = "bam", upper = 1))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = NA, upper = 1))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = NULL, upper = 1))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = 0, upper = "bam"))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = 0, upper = NA))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = 0, upper = NULL))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = 1, upper = -1))
  expect_error(makeIntegerVectorParam(id = "x", len = 2, lower = c(1, 1), upper = c(0,0)))
})


test_that("discrete param", {
  f = function(x) 2 * x
  p = makeDiscreteParam(id = "x", values = list(a = "char", b = 2L, c = 2.2, d = f, "e"))
  expect_equal("discrete", p$type)
  expect_true(isFeasible(p, "char"))
  expect_true(isFeasible(p, 2L))
  expect_true(isFeasible(p, 2.2))
  expect_true(isFeasible(p, f))
  expect_true(isFeasible(p, "char"))
  expect_true(isFeasible(p, 2L))
  expect_true(isFeasible(p, 2))
  expect_true(isFeasible(p, 2.2))
  expect_true(isFeasible(p, function(x) 2 * x))
  expect_true(isFeasible(p, "e"))

  expect_true(!isFeasible(p, "a"))
  expect_true(!isFeasible(p, sum))
  expect_true(!isFeasible(p, NULL))

  expect_equal(p$lower, NULL)
  expect_equal(p$upper, NULL)

  ## Error conditions:
  expect_error(makeDiscreteParam(id = "x", values = list(a = 1, "a")), "names could not be guessed")
  expect_error(makeDiscreteParam(id = "x", values = list()), "No possible value")
})

test_that("discrete vec param", {
  f = function(x) 2 * x
  p = makeDiscreteVectorParam(id = "x", len = 2, values = list(a = "char", b = 2L, c = 2.2, d = f, "e"))
  expect_equal("discretevector", p$type)
  expect_true(isFeasible(p, list("char", "char")))
  expect_true(isFeasible(p, list("char", f)))
  expect_true(isFeasible(p, list(2.2, 2L)))

  expect_true(!isFeasible(p, list()))
  expect_true(!isFeasible(p, list("char")))
  expect_true(!isFeasible(p, list("a", "b")))
  expect_true(!isFeasible(p, list(NULL, NULL)))

  expect_equal(p$lower, NULL)
  expect_equal(p$upper, NULL)

  ## Error conditions:
  expect_error(makeDiscreteVectorParam(id = "x", len = 2, values = list(a = 1, "a")), "names could not be guessed")
  expect_error(makeDiscreteVectorParam(id = "x", len = 2, values = list()), "No possible value")
})



test_that("logic param", {
  p = makeLogicalParam(id = "x")
  expect_equal("logical", p$type)
  expect_true(isFeasible(p, TRUE))
  expect_true(isFeasible(p, FALSE))

  expect_true(!isFeasible(p, "bam"))
  expect_true(!isFeasible(p, 1L))
  expect_true(!isFeasible(p, 1))
  expect_true(!isFeasible(p, NULL))
})

test_that("logic vec param", {
  p = makeLogicalVectorParam(id = "x", len = 2)
  expect_equal("logicalvector", p$type)

  expect_true(isFeasible(p, c(TRUE, FALSE)))

  expect_true(!isFeasible(p, "bam"))
  expect_true(!isFeasible(p, TRUE))
  expect_true(!isFeasible(p, FALSE))
  expect_true(!isFeasible(p, NULL))
})

test_that("character param", {
  p = makeCharacterParam(id = "s")
  expect_equal("character", p$type)
  expect_true(isFeasible(p, collapse(sample(letters, 5L))))

  expect_true(!isFeasible(p, 1L))
  expect_true(!isFeasible(p, 1))
  expect_true(!isFeasible(p, NULL))
  expect_true(!isFeasible(p, factor("bam")))
})

test_that("character vec param", {
  p = makeCharacterVectorParam(id = "x", len = 2)
  expect_equal("charactervector", p$type)
  expect_true(isFeasible(p, c("a", "b")))
  expect_false(isFeasible(p, c(1, 1)))
  expect_false(isFeasible(p, "a"))
  expect_false(isFeasible(p, 1))

  p = makeCharacterVectorParam(id = "x", len = 2, cnames = c("x1", "x2"))
  expect_equal("charactervector", p$type)
  expect_true(isFeasible(p, c(x1 = "a", x2 = "b")))
  expect_false(isFeasible(p, c("a", "b")))
})

test_that("function param", {
  p = makeFunctionParam(id = "x")
  expect_equal("function", p$type)
  expect_true(isFeasible(p, identity))

  expect_true(!isFeasible(p, "bam"))
  expect_true(!isFeasible(p, 1L))
  expect_true(!isFeasible(p, 1))
  expect_true(!isFeasible(p, NULL))
})

test_that("untyped param", {
  p = makeUntypedParam(id = "x")
  expect_equal("untyped", p$type)
  expect_true(isFeasible(p, 1))
  expect_true(isFeasible(p, identity))
  expect_true(isFeasible(p, NULL))
})

test_that("param print works", {
  p = makeNumericParam(id = "x", lower = -1L, upper = 1)
  expect_output(print(p), "numeric")
})



