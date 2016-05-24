test_that("CalcICV throws error",
{
  expect_that(CalcICV(array(1:100, c(10, 10))), 
              throws_error("covariance matrix must be symmetric."))
  expect_that(CalcICV(RandomMatrix(10)), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used for ICV."))
})
test_that("CalcICV returns sensible results",
{
  data(iris)
  cor.matrix = cor(iris[,1:4])
  cov.matrix = cov(iris[,1:4])
  eVals <- eigen(cov.matrix, only.values = TRUE)$values
  ident = diag(rep(2, 10))
  expect_that(CalcICV(ident), equals(0))
  expect_that(CalcICV(cov.matrix), equals(sd(eVals)/mean(eVals)))
}
)