context("HydeUtilities")

data(PE, package="HydeNet")
Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE)

test_that("decisionOptions - dbern",
{
  expect_equal(decisionOptions("treat", Net),
               0:1)
})

test_that("decisionOptions - dcat",
{
  expect_equal(decisionOptions("death", Net),
               1:2)
})

test_that("validateParameters",
{
  expect_equal(validateParameters(list(lambda = 5), dist = "dpois"),
               c("lambda > 0" = TRUE))
})

test_that("validateParameters - use fromData() and fromFormula()",
{
  expect_equal(validateParameters(list(mu = fromData(), 
                                       tau = fromFormula()), 
                                  dist = "dnorm"),
               c("is.numeric(mu)" = TRUE, "tau >= 0" = TRUE))
})