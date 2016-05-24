# Unit tests on the DidacticBoost package
if(require(rpart) & require(testthat)){
  context("Testing the package")
  test_that("Testing functions", {
    k <- kyphosis
    k$Kyphosis <- factor(ifelse(k$Kyphosis == "present", 1L, -1L))
    expect_that(fit <- fitBoosted(Kyphosis ~ Age + Number + Start, data = k, iterations = 10), not(throws_error()))
    expect_true(is.boosted(fit))
    expect_that(predict(fit), not(throws_error()))
    expect_that(predict(fit, newdata = k), not(throws_error()))
  })
}
