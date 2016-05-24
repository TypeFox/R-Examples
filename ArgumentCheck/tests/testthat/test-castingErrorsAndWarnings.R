context("finishArgCheck Errors and Warnings")

test_that("Cast Appropriate Error",
{
  Check <- newArgCheck()
  addError("New Error",
           Check)
  expect_error(finishArgCheck(Check),
               "1: New Error")

})

test_that("Cast Appropriate Warning",
{
  Check <- newArgCheck()
  addWarning("New Warning",
             Check)
  expect_warning(finishArgCheck(Check),
                 "1: New Warning")
})

test_that("Cast Simultaneous Error and Warning",
{
  Check <- newArgCheck()
  addError("New Error",
           Check)
  addWarning("New Warning",
             Check)
  expect_warning(expect_error(finishArgCheck(Check),
                              "1: New Error"),
                 "1: New Warning")
})