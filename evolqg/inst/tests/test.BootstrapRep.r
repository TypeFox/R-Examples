test_that("BootstrapRep returns sensible results",
          {
              ind.data <- iris[,1:4]
              expect_that(BootstrapRep(ind.data, RandomSkewers, iterations = 50) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, MatrixCor, iterations = 50, correlation = T) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, KrzCor, iterations = 50) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, KrzCor, iterations = 50, correlation = T) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, RandomSkewers, iterations = 50) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, RandomSkewers, iterations = 50, correlation = T) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, MatrixCor, iterations = 50, correlation  = T) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, KrzCor, iterations = 50) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, KrzCor, iterations = 50, correlation = T) >= 0.98 , is_true())
              expect_that(BootstrapRep(cov(ind.data), RandomSkewers, iterations = 50), throws_error("input appears to be a matrix, use residuals."))
              expect_that(BootstrapRep(cor(ind.data), RandomSkewers, iterations = 50), throws_error("input appears to be a matrix, use residuals."))
          }
)
test_that("BootstrapR2 returns sensible results",
{
  ind.data <- iris[,1:4]
  expect_that(length(BootstrapR2(ind.data)), equals(1000))
  expect_that(length(BootstrapR2(ind.data, 10)), equals(10))
  results <- BootstrapR2(ind.data, 10)
  corrs = sapply(results, function(x) isTRUE(x < 1 & x > 0))
  expect_that(all(corrs), is_true())
  expect_that(BootstrapR2(RandomMatrix(10), 10), throws_error("input appears to be a matrix, use residuals."))
})

