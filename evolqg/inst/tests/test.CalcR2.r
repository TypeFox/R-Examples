test_that("CalcR2 returns correct value for correlation matrices",
          {
            data(iris)
            cor.matrix  = cor(iris[,1:4])
            cov.matrix  = cov(iris[,1:4])
            ident = diag(rep(1, 10))
            expect_that(CalcR2(ident), equals(0))
            expect_that(CalcR2(cor.matrix), equals(mean((cor.matrix[lower.tri(cor.matrix)])^2)))
            expect_that(CalcR2(cov.matrix), equals(mean((cor.matrix[lower.tri(cor.matrix)])^2)))
          }
)
