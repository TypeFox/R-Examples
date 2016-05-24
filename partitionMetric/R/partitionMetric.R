partitionMetric <-
function (B, C, beta=2.0) {

  ## check input args
  if (length (B) != length (C)) {
    stop ('length of B must equal length of C')
  }

  if (beta <= 1.0) {
    stop ('beta must be > 1.0')
  }

  Bfac <- factor(B)
  Cfac <- factor(C)
  
  Bsum <- sum(table(Bfac)^beta)
  Csum <- sum(table(Cfac)^beta)
  intSum <- sum(table(Bfac, Cfac)^beta)

  (Bsum + Csum - 2*intSum)/((1-2^(1-beta))*(length(B)^beta))

}

