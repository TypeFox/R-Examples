context("ClippedFT")

test_that("clippedFT works as expected",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      Y2 <- rank(Y1)/64
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25,0.5)
      
      cFT1 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE)
      cFT2 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=FALSE)
      
      # Compute cFT 'by hand':
      res1 <- matrix(0, nrow=64, ncol=2)
      res2 <- matrix(0, nrow=64, ncol=2)
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          ii <- complex(real = 0, imaginary = 1)
          res1[f,l] <- sum((Y2 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
          res2[f,l] <- sum((Y1 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
        }
      }
      
      V1 <- getValues(cFT1, frequencies=freq, levels=levels)
      V2 <- getValues(cFT2, frequencies=freq, levels=levels)
      
      expect_that(dim(V1),equals(c(64,2,1)))
      expect_that(dim(V2),equals(c(64,2,1)))
      expect_that(V1[,,1],equals(res1))
      expect_that(V2[,,1],equals(res2))
      
    }
)
