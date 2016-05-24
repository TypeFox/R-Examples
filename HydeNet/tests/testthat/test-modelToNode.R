context("modelToNode")

test_that("modelToNode: coxph should return an error",
{
  library(survival)
  test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
                stop=c(2,3,6,7,8,9,9,9,14,17), 
                event=c(1,1,1,1,1,1,1,0,0,0), 
                x=c(1,0,0,1,0,1,1,1,0,0)) 
  fit <- coxph(Surv(start, stop, event) ~ x, test2)
  expect_error(modelToNode(fit))
})

test_that("modelToNode: multinom",
{
  fit.gear <- multinom(gear ~ mpg + factor(am), data=mtcars)
  expect_that(modelToNode(fit.gear), not(throws_error()))
})
  