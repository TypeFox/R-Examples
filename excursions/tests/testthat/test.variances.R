context("Variances")
test_that("Variances", {
  data <- integration.testdata1()
  vars <- excursions.variances(data$L, max.threads=1)
  v = diag(solve(data$Q))
  expect_equal(vars,v,tolerance=1e-7)
})


test_that("Variances Q and L", {
  data <- integration.testdata1()
  v1 <- excursions.variances(L = data$L, max.threads=1)
  v2 <- excursions.variances(Q = data$Q, max.threads=1)
  expect_equal(v1,v2,tolerance=1e-7)
})
