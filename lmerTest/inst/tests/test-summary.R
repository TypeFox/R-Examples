library("testthat")
library("lmerTest")

context("lmer models and methods")
test_that("summary fails for ProblemSet. Hessian is not positive definite", {
  load(system.file("testdata","ProblemSet.RData",package="lmerTest"))
  m <- lmer(ExpScore~ChLevel+(1|Subject), data=ProblemSet)
  #summary(m)
  expect_is(m,"merModLmerTest")
  #expect_error(summary(m)," length of 'dimnames' [2] not equal to array extent") 
});
