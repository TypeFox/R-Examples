context("ClippedCov")

test_that("clippedCov works as expected",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      Y2 <- rank(Y1)/64
      n <- length(Y1)

      levels <- c(0.25,0.5)
      K <- length(levels)
      
      cCov <- clippedCov(Y1, levels.1 = levels)
      
      # Compute cCov 'by hand':
      res1 <- array(0, dim = c(n, K, K, 1))
      res2 <- array(0, dim = c(n, K, K, 1))
      for (k1 in 1:K) {
        X1 <- (Y2 <= levels[k1]) - levels[k1]
        for (k2 in 1:K) {
          X2 <- (Y2 <= levels[k2]) - levels[k2]
          for (k in 0:(n-1)) {
            res1[k+1,k1,k2,1] <- (1/n) * sum( ( X1[1:(n-k)] ) * ( X2[(1+k):n]  ) )
          }
        }
      }
      
      V1 <- getValues(cCov, levels.1 = levels)
      
      expect_that(dim(V1),equals(c(64,2,2,1)))
      expect_that(V1,equals(res1))
      
    }
)
