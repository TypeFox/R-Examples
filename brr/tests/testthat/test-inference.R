context("Inference")

test_that("Credibility intervals", {
  a <- 0.5; b <- 0; c <- 0.5; d <- 0; S <- 10; T <- 5; x <- 10; y <- 10
  intervals <- brr_intervals(x, y, S, T, a, b, c, d, intervals=c("left","equi-tailed"))
  expect_is(intervals, "list")
  expect_identical(unname(intervals$equi[1]), qpost_phi(.025, a, b, c, d, S, T, x, y))
  model <- Brr(a=a,b=b,c=c,d=d,S=S,T=T,x=x,y=y)
  ci <- confint(model, intervals=c("left","equi-tailed"))
  attributes(ci) <- NULL
  expect_identical(unname(intervals), ci)
})