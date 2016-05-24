context("id special chars")


test_that("problem / algo ids", {
  reg = makeTestRegistry()
  expect_error(addProblem(reg, id="."), "comply with")
  expect_error(addProblem(reg, id=".foo"), "comply with")
  expect_error(addProblem(reg, id="-foor"), "comply with")

  expect_error(addAlgorithm(reg, id="."), "comply with")
  expect_error(addAlgorithm(reg, id=".foo"), "comply with")
  expect_error(addAlgorithm(reg, id="-foobar"), "comply with")
})
