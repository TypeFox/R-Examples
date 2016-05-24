
test_that("fitprobit() returns the right kind of object", {
  toxdat <- data.frame(
    dose=c(0.05, 0.0625, 0.125, 0.25, 0.5, 1),
    ntot=rep(8, 6),
    nfx = c(0, 1, 4, 4, 6, 8))
  fit <- fitprobit(toxdat)
  expect_that(fit, is_a("glm"))
})
