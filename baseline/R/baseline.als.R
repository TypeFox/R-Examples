baseline.als <- function(spectra, lambda=6, p=0.05, maxit = 20){
## Eilers baseline correction for Asymmetric Least Squares
## Migrated from MATLAB original by Kristian Hovde Liland
## $Id: baseline.als.R 170 2011-01-03 20:38:25Z bhm $
#
# INPUT:
# spectra - rows of spectra
# lambda  - 2nd derivative constraint
# p       - regression weight for positive residuals
#
# OUTPUT:
# baseline  - proposed baseline
# corrected - baseline corrected spectra
# wgts      - final regression weights

  mn       <- dim(spectra)
  baseline <- matrix(0,mn[1],mn[2])
  wgts     <- matrix(0,mn[1],mn[2])

  # Sparse empty matrix (m x m)
  empt     <- as.matrix.csr(0,mn[2],mn[2])

  # Diagonal sparse matrix (m x m)
  speye    <- empt
  diag(speye) <- 1
  D  <- diff(speye,differences=2)
  DD <- 10^lambda*t(D)%*%D

  # Iterate through spectra
  for(i in 1:mn[1]){
    w  <- rep.int(1, mn[2])
    y  <- spectra[i,]

    # Iterate restriction and weighting
    for(it in 1:maxit){
      W <- empt
      diag(W) <- w

      # Restricted regression
      z <- solve(W+DD,w*y)
      w_old <- w

      # Weights for next regression
      w <- p * (y > z) + (1 - p) * (y < z)
      sw <- sum(w_old != w)

      # Break if no change from last iteration
      if(sw == 0)
        break;
    }
    baseline[i,] <- z
    wgts[i,] <- w
  }
  corrected <- spectra-baseline
  list(baseline=baseline,corrected=corrected,wgts=wgts)
}
