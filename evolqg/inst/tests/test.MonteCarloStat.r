test_that("MonteCarloR2 returns sensible results",
          {
            expect_that(length(MonteCarloR2(diag(rep(2, 10)), 10)), equals(1000))
            expect_that(length(MonteCarloR2(diag(rep(2, 10)), 10, 10)), equals(10))
            results <- MonteCarloR2(diag(rep(2, 10)), 10)
            corrs = sapply(results, function(x) isTRUE(x < 1 & x > 0))
            expect_that(sum(corrs), equals(1000))
            expect_that(MonteCarloR2(RandomMatrix(10), 10), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling."))
          }
)

test_that("MonteCarloRep returns sensible results",
          {
            expect_that(MonteCarloRep(RandomMatrix(10), 30, RandomSkewers, iterations = 1), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling."))
            
            cov.matrix <- RandomMatrix(10, 1, 1, 10)
            expect_that(MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, iterations = 50) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, 30, MatrixCor, correlation = TRUE) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, 30, KrzCor) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, 30, KrzCor, correlation = TRUE) <= 1, is_true())
          }
)
test_that("MonteCarloStat throws error",
{
  expect_that(MonteCarloStat(array(1:100, c(10, 10)), 10, 10, RandomSkewers, cov), 
              throws_error("covariance matrix must be symmetric."))
})
test_that("MonteCarloRep returns correct results for scalled crappy matrices",
          {
           mat = RandomMatrix(35, 1, 1 ,10)
           crappy_mat = var(mvtnorm::rmvnorm(10, sigma = mat))
           rep_c = MonteCarloRep(crappy_mat, 10, PCAsimilarity)
           rep_sc = MonteCarloRep(100000*crappy_mat, 10, PCAsimilarity)
           expect_less_than(abs(rep_c - rep_sc), 0.05)
          }
)