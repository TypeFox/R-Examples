test_that("Autonomy returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              auto <- ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta)
              expect_that(Autonomy(cov.matrix, beta), equals(as.numeric(auto)))
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(Autonomy(singular.matrix), 
                             "matrix is singular, can't compute autonomy directly. Using nearPD, results could be wrong")
          }
)
test_that("ConditionalEvolvability returns correct restults",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              cond.evol <- (t (beta) %*% solve (cov.matrix) %*% beta)^(-1)
              expect_that(ConditionalEvolvability(cov.matrix, beta), equals(as.numeric(cond.evol)))
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(ConditionalEvolvability(singular.matrix), 
                             "matrix is singular, can't compute conditional evolvability directly. Using nearPD, results could be wrong")
          }
)
test_that("Constraints returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta = Normalize(rnorm(4))
              const = abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta))
              expect_that(Constraints(cov.matrix, beta), equals(as.numeric(const)))
          }
)
test_that("Evolvability returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              evol <- t (beta) %*% cov.matrix %*% beta
              expect_that(Evolvability(cov.matrix, beta), equals(as.numeric(evol)))
          }
)
test_that("Flexibility returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           flex <- t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta)
           expect_that(Flexibility(cov.matrix, beta), equals(as.numeric(flex)))
           expect_that(Flexibility(cov.matrix, beta) <=  1, is_true())
           expect_that(Flexibility(cov.matrix, beta) >= -1, is_true())
         }
)

test_that("Pc1Percent returns correct results",
          {
            cov.matrix = cov(matrix(rnorm(30*10), 30, 10))
            pc1 <- eigen(cov.matrix)$values[1]/sum(diag(cov.matrix))
            expect_that(Pc1Percent(cov.matrix), equals(pc1))
            expect_that(Pc1Percent(cov.matrix) <=  1, is_true())
            expect_that(Pc1Percent(cov.matrix) > 0, is_true())
          }
)
test_that("Respondability returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           repond <- Norm(cov.matrix %*% beta)
           expect_that(Respondability(cov.matrix, beta), equals(repond))
         }
)

test_that("MeanMatrixStatistics returns correct results",
          {
            set.seed(42)
            iris.stats <- MeanMatrixStatistics(cov(iris[,1:4]))
            test.values <- read.table("iris.stats")
            expect_that(iris.stats, is_equivalent_to(test.values[,1]))
          }
)
