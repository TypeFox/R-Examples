test_that("ExtendMatrix returns correct results",
          {
              cov.matrix = RandomMatrix(11)
              ret.dim = sample(2:11, 1)
              ext.matrix = ExtendMatrix(cov.matrix, ret.dim = ret.dim)
              expect_that(ext.matrix, is_a('list'))
              last.eval = eigen(cov.matrix)$values[ret.dim]
              eVals = eigen(ext.matrix[[1]])$values
              extended = sapply(eVals, function(x) isTRUE(all.equal(x, last.eval)))
              expect_that(sum(extended), equals(11-ret.dim+1))
              expect_that(sum(eVals>0), equals(length(eVals)))
              expect_that(ExtendMatrix(cov.matrix[1:9, 1:9]), gives_warning("matrix is too small"))
          }
)
