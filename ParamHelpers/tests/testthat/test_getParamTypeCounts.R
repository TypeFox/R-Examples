context("getParamTypeCounts")

test_that("getParamTypeCounts", {
  checkNonOccuringTypes = function(or, par.set) {
    sapply(setdiff(ParamHelpers:::getSupportedParamTypes(), getParamTypes(par.set)), function(type) {
      expect_equal(or[[type]], 0L)
    })
  }

  par.set = makeParamSet()
  checkNonOccuringTypes(getParamTypeCounts(par.set), par.set)

  par.set = makeParamSet(
    makeNumericParam("numeric1", lower = 0L, upper = 10L),
    makeIntegerParam("integer1", lower = 0L, upper = 5L),
    makeIntegerParam("integer2", lower = 0L, upper = 5L),
    makeDiscreteParam("discrete1", values = letters[1:5]),
    makeCharacterParam("character1")
  )

  or = getParamTypeCounts(par.set)
  expect_equal(or$numeric, 1L)
  expect_equal(or$integer, 2L)
  expect_equal(or$discrete, 1L)
  expect_equal(or$character, 1L)
  checkNonOccuringTypes(or, par.set)
})
