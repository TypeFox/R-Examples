context("multivariate")

test_that("multivariate functions work as expected",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      R1 <- rank(Y1)/64
      
      Y2 <- rnorm(64)
      R2 <- rank(Y2)/64
      
      Y <- cbind(Y1, Y2)
      R <- cbind(R1, R2)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25,0.5)
      
      cFT_R <- clippedFT(Y, frequencies = freq, levels = levels, isRankBased=TRUE)
      cFT_Y <- clippedFT(Y, frequencies = freq, levels = levels, isRankBased=FALSE)
      
      # Compute cFT 'by hand':
      res_Y <- array(0, dim = c(64, 2, 2))
      res_R <- array(0, dim = c(64, 2, 2))
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          for (i in 1:2) {
            ii <- complex(real = 0, imaginary = 1)
            res_Y[f,i,l] <- sum((Y[,i] <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
            res_R[f,i,l] <- sum((R[,i] <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))            
          }
        }
      }
      
      V_Y <- getValues(cFT_Y, frequencies=freq, levels=levels)
      V_R <- getValues(cFT_R, frequencies=freq, levels=levels)
      
      expect_that(dim(V_Y),equals(c(64,2,2,1)))
      expect_that(dim(V_R),equals(c(64,2,2,1)))
      expect_that(V_Y[,,,1],equals(res_Y))
      expect_that(V_R[,,,1],equals(res_R))
      
      # Compute qPG 'by hand':
      res_Y_qPG <- array(0, dim = c(64, 2, 2, 2, 2))
      res_R_qPG <- array(0, dim = c(64, 2, 2, 2, 2))
      
      for (l1 in 1:length(levels)) {
        for (l2 in 1:length(levels)) {
          for (i1 in 1:2) {
            for (i2 in 1:2) {
              res_Y_qPG[,i1,l1,i2,l2] <- res_Y[,i1,l1] * Conj(res_Y[,i2,l2]) / (2*pi*64)
              res_R_qPG[,i1,l1,i2,l2] <- res_R[,i1,l1] * Conj(res_R[,i2,l2]) / (2*pi*64)
            }
          }
        }
      }     
  
      qPG_R <- quantilePG(Y, frequencies = freq, levels.1 = levels, isRankBased=TRUE, type="clipped")
      qPG_Y <- quantilePG(Y, frequencies = freq, levels.1 = levels, isRankBased=FALSE, type="clipped")
      
      V_R <- getValues(qPG_R)
      V_Y <- getValues(qPG_Y)
      
      expect_that(dim(V_Y),equals(c(64,2,2,2,2,1)))
      expect_that(dim(V_R),equals(c(64,2,2,2,2,1)))
      expect_that(V_Y[,,,,,1],equals(res_Y_qPG))
      expect_that(V_R[,,,,,1],equals(res_R_qPG))
      
      ###################################
      # perform smoothing 'by hand':
      ###################################
      
      res_Y_sPG <- array(0, dim = c(64, 2, 2, 2, 2))
      res_R_sPG <- array(0, dim = c(64, 2, 2, 2, 2))
      
      weight <- kernelWeight(N=64, W=W0, bw=0.1)
      
      #sPG <- smoothedPG(Y, frequencies = freq, levels.1 = levels, weight=weight)
      sPG_R <- smoothedPG(qPG_R, weight=weight, frequencies=freq)
      sPG_Y <- smoothedPG(qPG_Y, weight=weight, frequencies=freq)
      
      # perform smoothing 'by hand':
      
      # weight function
      Wn <- function(x) {
        start <- ceiling(-weight@bw-x/(2*pi))
        end <- floor(weight@bw-x/(2*pi))
        result <- 0
        for (j in start:end) {
          result <- result + weight@bw^(-1) * weight@W(weight@bw^(-1)*(x+2*pi*j))
        }
        result
      }
      
      ref_R <- array(0, dim=c(64,2,2,2,2))
      ref_Y <- array(0, dim=c(64,2,2,2,2))
      for (f in 1:length(freq)) {
        for (s in 1:63) {
          ref_R[f,,,,] <- ref_R[f,,,,] + Wn(freq[f]-2*pi*s/64) * res_R_qPG[s+1,,,,]
          ref_Y[f,,,,] <- ref_Y[f,,,,] + Wn(freq[f]-2*pi*s/64) * res_Y_qPG[s+1,,,,]
        }
        ref_R[f,,,,] <- 2*pi*ref_R[f,,,,]/ (64 * (weight@env$Wnj[c(64,1:63)])[f])
        ref_Y[f,,,,] <- 2*pi*ref_Y[f,,,,]/ (64 * (weight@env$Wnj[c(64,1:63)])[f])
      }
      
      
      S_R <- getValues(sPG_R)
      S_Y <- getValues(sPG_Y)
      
      expect_that(dim(S_R),equals(c(64,2,2,2,2,1)))
      expect_that(S_R[,,,,,1],equals(ref_R))
      expect_that(dim(S_Y),equals(c(64,2,2,2,2,1)))
      expect_that(S_Y[,,,,,1],equals(ref_Y))
      
    }
)


test_that("multivariate functions work as expected (check against ref-results.rdata",{

      source("load-ref.R")
      
      set.seed(2581)
      Y1 <- ts1(16)
      Y2 <- ts1(16)
      
      Y <- matrix(c(Y1, Y2), ncol=2)
      
      lev.ok.all <- c(0.25,0.5,0.75)
      weight = kernelWeight(W=W0, N=64, bw=0.2)
      
      sPG.fft <- smoothedPG(Y, levels.1=lev.ok.all, weight = weight, type="clipped")
      
      W.fft.ci.mult <- getPointwiseCIs(sPG.fft)
      W.fft.coh.ci.mult <- getPointwiseCIs(sPG.fft, quantity = "coherency")
      W.fft.cohsq.ci.mult <- getPointwiseCIs(sPG.fft, quantity = "coherence")
      
      W.fft.sd.mult <- getSdNaive(sPG.fft)
      W.fft.coh.sd.1.mult <- getCoherencySdNaive(sPG.fft, type="1")
      W.fft.coh.sd.2.mult <- getCoherencySdNaive(sPG.fft, type="2")
      
      
      expect_equal(W.fft.ci.mult, W.fft.ci.mult.ref)
      expect_equal(W.fft.coh.ci.mult.ref, W.fft.coh.ci.mult)
      expect_equal(W.fft.cohsq.ci.mult.ref, W.fft.cohsq.ci.mult)
      expect_equal(W.fft.sd.mult.ref, W.fft.sd.mult.ref)
      expect_equal(W.fft.coh.sd.1.mult.ref, W.fft.coh.sd.1.mult.ref)
      expect_equal(W.fft.coh.sd.2.mult.ref, W.fft.coh.sd.2.mult.ref)
      
    }
)
