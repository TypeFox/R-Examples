context("Test summary.ltmle") 

test_that("OR and RR are not calculated when Y is not binary", {
	n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=list(1, 0), IC.variance.only=TRUE)

  s <- summary(r1)

  expect_that(names(s$effect.measures), equals(c("treatment", "control", "ATE")))

  t1 <- ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=list(1, 0), survivalOutcome=FALSE, estimate.time = FALSE)
  u <- summary(t1)

  expect_that(names(u$effect.measures), equals(c("treatment", "control", "ATE", "RR", "OR")))

})

test_that("CI is truncated if outcome is binary", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(W * 10)
  Y <- rbinom(n, 1, 0.99)
  data <- data.frame(W, A, Y)
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, IC.variance.only = FALSE, estimate.time = FALSE)
  s1 <- summary(r1)
  expect_less_than(max(s1$treatment$CI), 1.0001)
})

test_that("CI is not truncated if outcome is not binary", {
  n <- 100
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rnorm(n, mean = 5))
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, IC.variance.only = TRUE, estimate.time = FALSE)
  s1 <- summary(r1)
  expect_more_than(max(s1$treatment$CI), 1)
})