###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

# method of Theorem 1.1 in Bogdan et. al. (2015)
lambdaBH <- function(fdr, n.group) {
  lambda.BH <- rep(NA,n.group)
  for (i in 1:n.group){
    lambda.BH[i] <- qnorm(1-(i*fdr)/(2*n.group))
  }

  return(lambda.BH)
}

# method of Section 3.2.2 in Bogdan et. al. (2015)
lambdaGaussian <- function(fdr, n.group, n.obs) {
  lambda.BH <- lambdaBH(fdr=fdr, n.group=n.group)

  omegafun <- function(k) { return(1/(n.obs-k-1)) }

  lambda.G <- rep(NA,n.group)
  lambda.G[1] <- lambda.BH[1]
  for (i in 2:min(n.group, n.obs-2)) {
    lambda.G[i] <- lambda.BH[i] * sqrt( 1 + omegafun(i-1) * sum(lambda.G[1:(i-1)]^2) )
  }

  lambda.G.min <- min(lambda.G, na.rm=TRUE)
  min.ind <- which(lambda.G==lambda.G.min)
  lambda.G[min.ind:n.group] <- lambda.G.min

  return(lambda.G)
}

# method introduced in Gossmann et. al. (2015)
lambdaGaussianMC <- function(fdr, n.group, group.id, A, n.MC, MC.reps) {
  lambda.BH <- lambdaBH(fdr=fdr, n.group=n.group)

  mA <- matrix(NA, nrow(A), n.group)
  for (i in 1:n.group) {
    mA[ , i] <- apply(as.matrix(A[ , group.id[[i]] ]), 1, mean)
    mA[ , i] <- mA[ , i] / sqrt(sum(mA[ , i]^2))
  }

  # Monte Carlo corrections for lambda.BH
  lambda.MC <- lambdaMC(lambda.BH, mA, n.MC, MC.reps)
  lambda.MC <- c(lambda.MC, rep(lambda.MC[n.MC], n.group-n.MC))

  return(lambda.MC)
}

# lambdas of Theorem 2.5 and equation (2.14) in Brzyski et. al. (2015)
lambdaChiOrtho <- function(fdr, n.group, group.sizes, wt, method) {
  lambda.max <- rep(NA, n.group)
  lambda.min <- rep(NA, n.group)

  for (i in 1:n.group) {
    qchi.seq <- rep(NA, n.group)
    for (j in 1:n.group) {
      qchi.seq[j] <- sqrt(qchisq(1 - fdr*i/n.group, df=group.sizes[j])) / wt[j]
    }
    lambda.max[i] <- max(qchi.seq)
    lambda.min[i] <- min(qchi.seq)
  }

  # stop here if method is "chiOrthoMax"
  if (method=="chiOrthoMax") return(lambda.max)

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  lambda.mean <- rep(NA, n.group)
  for (k in 1:n.group) {
    if (lambda.min[k] == lambda.max[k]) {
      lambda.mean[k] <- lambda.max[k]
    } else {
      # compute inverse of cdfMean at 1-fdr*k/n.group
      cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*k/n.group)),
                             lower = lambda.min[k], upper = lambda.max[k], extendInt="yes")
      lambda.mean[k] <- cdfMean.inv$root
    }
  }

  return(lambda.mean)
}

# Procedure 1 in Brzyski et. al. (2015)
lambdaChiEqual <- function(fdr, n.obs, n.group, m, w) {
  lambda.chi    <- rep(NA, n.group)
  lambda.chi[1] <- sqrt(qchisq(1 - fdr / n.group, df=m)) / w
  # make sure that the denominator in s is non-zero
  imax <- floor((n.obs - 1) / m + 1)
  imax <- min(n.group, imax)
  for (i in 2:imax) {
    s <- (n.obs - m*(i-1)) / n.obs + (w^2 * sum(lambda.chi[1:(i-1)]^2)) / (n.obs - m*(i-1) - 1)
    s <- sqrt(s)
    lambda.tmp <- (s/w) * sqrt(qchisq(1 - fdr * i / n.group, df=m))
    if (lambda.tmp <= lambda.chi[i-1]) {
      lambda.chi[i] <- lambda.tmp
    } else {
      for (j in i:imax) { lambda.chi[j] <- lambda.chi[i-1] }
      break
    }
  }
  lambda.chi[imax:n.group] <- lambda.chi[imax]

  return(lambda.chi)
}

# Procedure 2 in Brzyski et. al. (2015)
lambdaChiMean <- function(fdr, n.obs, n.group, group.sizes, wt) {
  lambda.chi.mean <- rep(NA, n.group)

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  # get upper and lower bounds for lambda.chi.mean[1]
  qchi.seq <- rep(NA, n.group)
  for (j in 1:n.group) {
    qchi.seq[j] <- sqrt(qchisq(1 - fdr/n.group, df=group.sizes[j])) / wt[j]
  }
  upperchi <- max(qchi.seq)
  lowerchi <- min(qchi.seq)

  if (upperchi == lowerchi) {
    lambda.chi.mean[1] <- upperchi
  } else {
    lambda.chi.mean[1] <- uniroot(function(y) (cdfMean(y) - (1-fdr/n.group)),
                                  lower = lowerchi, upper = upperchi, extendInt="yes")$root
  }

  # get lambda.chi.mean[2:n.group]
  for (i in 2:n.group) {
    s <- rep(NA, n.group)
    for (j in 1:n.group) {
      # prevent division by 0
      if ((n.obs - group.sizes[j]*(i-1) - 1) <= 0) {
        for (k in i:n.group) { lambda.chi.mean[k] <- lambda.chi.mean[i-1] }
        break
      }

      s[j] <- (n.obs - group.sizes[j]*(i-1)) / n.obs + 
        (wt[j]^2 * sum(lambda.chi.mean[1:(i-1)]^2)) / (n.obs - group.sizes[j]*(i-1) - 1)
      s[j] <- sqrt(s[j])
    }

    cdfMean <- function(x) {
      pchi.seq <- rep(NA, n.group)
      for (j in 1:n.group) {
        # Procedure 2 in Brzyski et. al. (2015) has a typo at this point
        pchi.seq[j] <- pchisq((wt[j]/s[j] * x)^2, df=group.sizes[j])
      }
      return(mean(pchi.seq))
    }

    cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*i/n.group)),
                           lower = 0, upper = upperchi, extendInt="upX")$root

    if (cdfMean.inv <= lambda.chi.mean[i-1]) {
      lambda.chi.mean[i] <- cdfMean.inv 
    } else {
      for (j in i:n.group) { lambda.chi.mean[j] <- lambda.chi.mean[i-1] }
      break
    }
  }

  return(lambda.chi.mean)
}

# approximates the variance of (2.25) in Brzyski et. al. (2015) via Monte Carlo
lambdaChiMC <- function(fdr, X, y, group.id, wt, n.MC, MC.reps) {
  n.group     <- length(group.id)
  group.sizes <- sapply(group.id, FUN=length)
  lambda.MC   <- vector()

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  # get upper and lower bounds for lambda.MC[1]
  qchi.seq <- rep(NA, n.group)
  for (j in 1:n.group) {
    qchi.seq[j] <- sqrt(qchisq(1 - fdr/n.group, df=group.sizes[j])) / wt[j]
  }
  upperchi <- max(qchi.seq)
  lowerchi <- min(qchi.seq)

  if (upperchi == lowerchi) {
    lambda.MC[1] <- upperchi
  } else {
    lambda.MC[1] <- uniroot(function(y) (cdfMean(y) - (1-fdr/n.group)),
                            lower = lowerchi, upper = upperchi, extendInt="yes")$root
  }

  # get lambda.MC[2:n.MC]
  for (i in 2:n.MC) {
    s <- lambdaChiMCAdjustment(y=y, X=X, group_id=group.id, lambda=lambda.MC,
                               w=wt, number_of_drawings=MC.reps)

    cdfMean <- function(x) {
      pchi.seq <- rep(NA, n.group)
      for (j in 1:n.group) {
        pchi.seq[j] <- pchisq((wt[j]/s * x)^2, df=group.sizes[j])
      }
      return(mean(pchi.seq))
    }

    cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*i/n.group)),
                           lower = 0, upper = upperchi, extendInt="upX")$root

    if (cdfMean.inv <= lambda.MC[i-1]) {
      lambda.MC[i] <- cdfMean.inv 
    } else {
      lambda.MC[i] <- lambda.MC[i-1] 
    }
  }

  # get lambda.MC[n.MC:n.group]
  lambda.MC[(n.MC+1):n.group] <- lambda.MC[n.MC]
  
  return(lambda.MC)
}
