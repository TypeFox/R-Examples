context("repairPoint")

test_that("repairPoint", {

  ps = makeParamSet(
    makeNumericParam("num1", lower = 1, upper = 5),
    makeIntegerParam("num2", lower = 1, upper = 10),
    makeNumericVectorParam("num3", len = 2L, lower = 1, upper = 5),
    makeDiscreteParam("disc1", values = letters[1:10])
  )

  # sample values violating bounds for numeric params
  x = list(num1 = -2, num2 = 11, num3 = c(2, 10), disc1 = "b")

  x.repaired = repairPoint(ps, x)

  expect_true(x.repaired$num1 == 1)
  expect_true(x.repaired$num2 == 10)
  expect_true(all(x.repaired$num3 == c(2, 5)))
  expect_true(x.repaired$disc1 == x$disc1)
})


