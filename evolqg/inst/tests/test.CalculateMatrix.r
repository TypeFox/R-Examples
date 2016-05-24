test_that("CalculateMatrix returns correct result",
          {
            data(iris)
            options(contrasts=c("contr.sum","contr.poly"))
            iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
            cov.matrix = var(iris.lm$residuals)*((dim(iris.lm$residuals)[1]-1)/iris.lm$df.residual)
            expect_that(CalculateMatrix(iris.lm), equals(cov.matrix))
            expect_that(CalculateMatrix(lm(as.matrix(iris[,1:4])~1)), equals(cov(iris[,1:4])))
          }
)
