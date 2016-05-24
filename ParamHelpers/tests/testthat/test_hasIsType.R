context("is{Discrete, Integer, Numeric, Character} and has{Discrete, Integer, Numeric, Character}")

test_that("is{Discrete, Integer, Numeric, Character} and has{Discrete, Integer, Numeric, Character}", {
  par.set.empty = makeParamSet()
  has.methods = c(hasInteger, hasDiscrete, hasNumeric, hasCharacter)
  is.methods = c(isInteger, isDiscrete, isNumeric, isCharacter)

  lapply(has.methods, function(fun) expect_false(fun(par.set.empty)))
  lapply(is.methods, function(fun) expect_true(fun(par.set.empty)))

  par.set.mixed = makeParamSet(
    makeNumericParam("u", lower = 1),
    makeIntegerParam("v", lower = 1, upper = 2),
    makeDiscreteParam("w", values = 1:2),
    makeLogicalParam("x"),
    makeDiscreteVectorParam("y", len = 2, values = c("a", "b"))
  )

  par.set.discrete = makeParamSet(
    makeDiscreteParam("discr1", values = 1:2),
    makeDiscreteParam("discr2", values = letters[1:3]),
    makeDiscreteVectorParam("y", len = 2, values = c("x", "y"))
  )

  par.set.character = makeParamSet(
    makeCharacterParam("char1"),
    makeCharacterVectorParam("chars1", len = 2)
  )

  expect_true(isDiscrete(par.set.discrete))
  expect_false(isInteger(par.set.discrete))
  expect_false(isNumeric(par.set.discrete))
  expect_false(hasNumeric(par.set.discrete))
  expect_false(hasInteger(par.set.discrete))
  expect_true(isCharacter(par.set.character))
  expect_false(hasInteger(par.set.character))

  expect_false(isDiscrete(par.set.mixed))
  expect_false(isInteger(par.set.mixed))
  expect_false(isNumeric(par.set.mixed))
  expect_false(isCharacter(par.set.mixed))
  expect_true(hasNumeric(par.set.mixed))
  expect_true(hasInteger(par.set.mixed))
  expect_false(hasCharacter(par.set.mixed))
})
