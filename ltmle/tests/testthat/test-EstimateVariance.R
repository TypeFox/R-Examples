context("Test variance estimation")

test_that("TMLE based variance estimate > IC based variance estimate under positivity",{
  niter <- 20
  std.dev.ratio <- numeric(niter)
  for (i in 1:niter) {
    n <- 1000
    W <- rnorm(n)
    A <- rbinom(n, 1, plogis(5*W))
    Y <- rbinom(n, 1, 0.5)
    r1 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
    r2 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, IC.variance.only=TRUE)
    std.dev.ratio[i] <- summary(r1)$treatment$std.dev / summary(r2)$treatment$std.dev
  }
  expect_true(mean(std.dev.ratio) > 1.2)
})