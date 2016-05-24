context("extractParamSetFromSooFunction")

# FIXME: test disabled because of annying soobench interface change

if (FALSE) {

test_that("extraction of soobench function parmameter set works", {
  library(soobench)
  soo.fun = ackley_function(4)
  par.set = extractParamSetFromSooFunction(soo.fun)
  expect_true(all(getUpper(par.set) == upper_bounds(soo.fun)))
  expect_true(all(getLower(par.set) == lower_bounds(soo.fun)))
  expect_equal(getParamTypeCounts(par.set)$numeric, number_of_parameters(soo.fun))
  expect_true(isNumeric(par.set))
})

}
