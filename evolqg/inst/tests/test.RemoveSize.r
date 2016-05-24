test_that("RemoveSize returns correct results",
          {
            cov.matrix <- cov(iris[,1:4])
            eigen.size <- eigen(cov.matrix)
            eigen.no.size <- eigen(RemoveSize(cov.matrix))
            dot.eigens <- abs(eigen.size$vectors[,1] %*% eigen.no.size$vectors[,4])[1]
            expect_that(dot.eigens, equals(1))
            expect_that(eigen.no.size$values[1:3], is_equivalent_to(eigen.size$values[2:4]))
            expect_that(eigen.no.size$values[4], equals(0))
          }
)
