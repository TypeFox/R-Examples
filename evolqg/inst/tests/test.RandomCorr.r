test_that("RandomCorr returns sensible results",
          {
            mat.list <- lapply(as.list(1:100), function(x) RandomMatrix(10))
            mean.corr <- mean(MatrixCor(mat.list)[upper.tri(array(0, c(100, 100)))])
            expect_that(dim(mat.list[[1]]), equals(c(10, 10)))
            expect_that(sum(diag(mat.list[[1]])), equals(10))
            expect_that(abs(mean.corr) < 0.01, is_true())
          }
)

test_that("RandomMatrix returns sensible results",
          {
            mat.list <- RandomMatrix(10, 10)
            expect_that(dim(mat.list[[1]]), equals(c(10, 10)))
            expect_that(length(mat.list), equals(10))
            expect_that(sum(diag(mat.list[[1]])), equals(10))
            rand.mat <- RandomMatrix(10, 1, 1, 2)
            expect_that(sum(diag(rand.mat))> 10, is_true())
            expect_that(sum(diag(rand.mat))< 20, is_true())
            mat.list.fixed.var <- RandomMatrix(3, 2, variance = c(1, 2, 3))
            expect_that(diag(mat.list.fixed.var[[1]]), equals(diag(mat.list.fixed.var[[2]])))
            expect_that(diag(mat.list.fixed.var[[1]]), equals(c(1, 2, 3)))
          }
)
