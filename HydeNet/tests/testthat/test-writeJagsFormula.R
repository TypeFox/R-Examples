context("writeJagsFormula")

test_that("writeJagsFormula: Poisson Regression",
{
  fit <- glm(gear ~ mpg + am, data = mtcars, family = poisson)
  expect_that(writeJagsFormula(fit, c("gear", "mpg", "am")),
              not(throws_error()))
})

test_that("writeJagsFormula: Multinomial Regression",
{
  fit.gear <- multinom(gear ~ mpg + factor(am), data=mtcars)
  expect_that(writeJagsFormula(fit.gear), not(throws_error()))
})
