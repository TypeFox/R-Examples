.mimm <- function(R, LblCut, NVC, NF, Iter) {

      NItems <- NVC

  # Compute item-factor and fac-fac covariances by calculating sums.
  #   First calculate the item-fac covariances in I <- 1..NItems by summing
  #   within-cluster columns. Once these are computed, calculate the
  #   factor-factor covariances in I <- NItems+1..NF part of first loop by
  #   applying the same formula to the recently computed item-factor
  #   covariances.

  # R[L,I): at this point these are covariances, not correlations.
  # L: ordinal position of the Jth fac in the list of items and factors.

      for (I in 1:(NItems+NF)) {
        for (J in 1:NF) {
          L <- NItems + J
          R[I,L] <- 0.0
          for (K in LblCut[J,1]:LblCut[J,2])
            R[I,L] <- R[I,L] + R[I,K]
          R[L,I] <- R[I,L]
        }
      }

  # Compute coefficient alpha for each group.
  # N: Number of items in the Jth group.
  # R[L,L): Jth factor variance, computed in previous paragraph.
  # DSum: sum of the diagonal of the Jth group.
  # RMean: mean of the off-diagonal correlations of Jth group.
      
      Alpha <- numeric(length=NF)

      for (J in 1:NF) {
        L <- NItems + J
        N <- LblCut[J,2] - LblCut[J,1] + 1
        if (N == 1)
          Alpha[J] <- 1.0
        else {
          DSum <- 0.0
          for (I in LblCut[J,1]:LblCut[J,2]) {
            DSum <- DSum + R[I,I]
          }
          XN <- as.numeric(N)
          RMean <- (R[L,L]-DSum) / (XN*(XN-1.0))
          Alpha[J] <- (XN*RMean) / ((XN-1)*RMean+1.0)
        }
      }

  # Compute communalities of the Jth group by iterating.
  # Iter: Number of iterations.
  # OldFF: Factor variance at the beginning of current iteration.
  # OldII: Communality at the beginning of current iteration.

      if (Iter >= 0) {

        for (J in 1:NF) {
          if (Alpha[J] > 0) {
            for (It in 1:Iter) {
              L <- NItems + J
              OldFF <- R[L,L]
              R[L,L] <- 0.0
              for (I in LblCut[J,1]:LblCut[J,2]) {
                OldII <- R[I,I]
                R[I,I] <- R[I,L]**2 / OldFF
                R[I,L] <- R[I,L] + (R[I,I] - OldII)
                R[L,L] <- R[L,L] + R[I,L]
              }
            }
          }
        }

      }

  # Standardize the item-factor and factor-factor covariances.
  # Since the items are already standardized, item-factor covariances are
  #   standardized by dividing by FacSD. This is accomplished in a loop
  #   as well as beginning the standardization of the factor-factor
  #   covariances.  The second division by the remaining FacSD for the
  #   factor-factor covariances is accomplished in a following loop.

      for (J in 1:NF) {
        L <- NItems + J
        FacSD <- sqrt(R[L,L])
        for (I in 1:L) {
          R[I,L] <- R[I,L] / FacSD
          R[L,I] <- R[I,L]
        }
        for (I in J:NF) {
          LI <- NItems + I
          R[L,LI] <- R[L,LI] / FacSD
          R[LI,L] <- R[L,LI]
        }
      }

  # Compute the Joreskog (1971) reliability coef, which is Omega.
  # Omega uses each items communality as a reliability estimate.
  # Alpha assumes equal item reliabilities.

      Omega <- numeric(length=NF)
      for (IFac in 1:NF) {
        FI <- NItems + IFac
        SumLam <- 0
        SumUnq <- 0
        for (J in LblCut[IFac,1]:LblCut[IFac,2]) {
          Lam <- R[FI,J]
          Unique <- 1 - Lam**2
          SumLam <- SumLam + Lam
          SumUnq <- SumUnq + Unique
        }
        Omega[IFac] <- (SumLam**2) / (SumLam**2 + SumUnq)
        if (Iter == 0) Omega[IFac] <- 1.0
      }

  return(list(R=R, Alpha=Alpha, Omega=Omega))

}

