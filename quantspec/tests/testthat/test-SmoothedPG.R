context("SmoothedPG")

test_that("smoothedPG works as expected- compare to naively computed estimator",{
      
      set.seed(2581)
      Y <- rnorm(64)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25,0.5)
      
      qPG <- quantilePG(Y, frequencies = freq, levels.1 = levels)
      
			weight <- kernelWeight(N=64, W=W0, bw=0.1)
			#weight <- kernelWeight(W=W0, bw=0.1)
			
			#sPG <- smoothedPG(Y, frequencies = freq, levels.1 = levels, weight=weight)
      sPG <- smoothedPG(qPG, weight=weight, frequencies=freq)
      
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
      
      Q <- getValues(qPG)

      ref <- array(0, dim=c(64,2,2,1))
      for (f in 1:length(freq)) {
        for (s in 1:63) {
          ref[f,,,] <- ref[f,,,] + Wn(freq[f]-2*pi*s/64) * Q[s+1,,,] 
        }
        ref[f,,,] <- 2*pi*ref[f,,,]/ (64 * (weight@env$Wnj[c(64,1:63)])[f])
      }
      

      S <- getValues(sPG)

      expect_that(dim(S),equals(c(64,2,2,1)))
      expect_that(S,equals(ref))
    }
)


test_that("smoothedPG (with SpecDistrWeight) works as expected- compare to naively computed estimator",{
			
			set.seed(2581)
			
			Y <- rnorm(64)
			
			freq <- 2*pi*(0:63)/64
			levels <- c(0.25,0.5)
			
			qPG <- quantilePG(Y, frequencies = freq, levels.1 = levels)
			weight <- specDistrWeight()
			sPG <- smoothedPG(qPG, weight=weight, frequencies=freq)
			
			# perform smoothing 'by hand':

			Q <- getValues(qPG)
			
			ref <- array(0, dim=c(64,2,2,1))
			for (f in 1:length(freq)) {
				for (s in 1:63) {
					ref[f,,,] <- ref[f,,,] + (freq[f] >= 2*pi*s/64) * Q[s+1,,,] 
				}
				ref[f,,,] <- 2*pi*ref[f,,,] / 64
			}
			
			
			S <- getValues(sPG)
			
			expect_that(dim(S),equals(c(64,2,2,1)))
			expect_that(S,equals(ref))
		}
)

test_that("smoothedPG works as expected for various levels",{
      
      source("load-ref.R")
      
      lev.ok.all <- c(0.25,0.5,0.75)
      lev.ok.1   <- c(0.25)
      lev.ok.2   <- c(0.5,0.75)
      lev.err.1  <- c("non numeric",0.5)
      lev.err.2  <- c(0.5,1.5)
      
      # Check whether it works for all levels:
      W.qr.ref.1 <- array(W.qr.ref[,,,1], dim=c(64,3,3,1))
      W.fft.ref.1 <- array(W.fft.ref[,,,1], dim=c(64,3,3,1))
      
      weight = kernelWeight(W=W0, N=64, bw=0.2)
			#weight = kernelWeight(W=W0, bw=0.2)
      
      sPG.qr <- smoothedPG(Y, levels.1=lev.ok.all, weight = weight, type="qr")
      W.qr <- getValues(sPG.qr)      
      expect_that(dim(W.qr),equals(c(64,3,3,1)))
      expect_that(W.qr,equals(W.qr.ref.1))
      
      sPG.fft <- smoothedPG(Y, levels.1=lev.ok.all, weight = weight, type="clipped")
      W.fft <- getValues(sPG.fft)
      expect_that(dim(W.fft),equals(c(64,3,3,1)))
      expect_that(W.fft,equals(W.fft.ref.1))

      W.fft.sd <- getSdNaive(sPG.fft)
      expect_that(dim(W.fft.sd),equals(c(64,3,3)))
      expect_that(W.fft.sd,equals(W.fft.sd.ref))
      
    })