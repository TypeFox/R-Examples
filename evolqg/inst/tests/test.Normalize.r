test_that("Normalize returns correct results",
          {
            beta <- rep(1, 10)
            expect_that(Normalize(beta), equals(rep(1/sqrt(10), 10)))
            expect_that(Norm(Normalize(beta)), equals(1))
          }
)

test_that("Norm returns correct results",
          {
            beta <- rep(1, 10)
            expect_that(Norm(beta), equals(sqrt(sum((rep(1, 10))^2))))
            expect_that(Norm(rep(1/sqrt(10), 10)), equals(1))
          }
)
