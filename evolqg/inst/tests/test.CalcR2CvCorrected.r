test_that("CalcR2CvCorrected returns sensible results on residuals",
          {
              results <- CalcR2CvCorrected(iris[,1:4], iterations = 100)
              expect_that(length(results), equals(3))
              expect_that(results, is_a('list'))
              expect_that(results[[1]], is_a('numeric'))
              expect_that(results[[2]], is_a('list'))
              expect_that(results[[3]], is_a('data.frame'))
              expect_that(names(results), equals(c("adjusted.integration.index",
                                                      "models", "dist")))
              expect_that(colnames(results[[3]]), equals(c("r2", "eVals_cv", "mean_cv")))
              expect_that(sum(sapply(results[[3]][,1], function(x) !(x > 0 & x < 1))), equals(0))
              expect_that(length(results[[1]]), equals(2))
              expect_that(length(results[[2]]), equals(2))
              expect_that(dim(results[[3]]), equals(c(100, 3)))
          })

test_that("CalcR2CvCorrected returns sensible results on models",
          {
              iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
              results <- CalcR2CvCorrected(iris.lm, iterations = 100)
              expect_that(length(results), equals(3))
              expect_that(results, is_a('list'))
              expect_that(results[[1]], is_a('numeric'))
              expect_that(results[[2]], is_a('list'))
              expect_that(results[[3]], is_a('matrix'))
              expect_that(names(results), equals(c("adjusted.integration.index",
                                                      "models", "dist")))
              expect_that(colnames(results[[3]]), equals(c("r2", "eVals_cv", "mean_cv")))
              expect_that(sum(sapply(results[[3]][,1], function(x) !(x > 0 & x < 1))), equals(0))
              expect_that(length(results[[1]]), equals(2))
              expect_that(length(results[[2]]), equals(2))
              expect_that(dim(results[[3]]), equals(c(100, 3)))
          })
