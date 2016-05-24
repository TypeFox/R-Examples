gsd.bdryconstant <- function(ifrac, eprob = 0.05, delta = 0.5, alternative = c("two.sided", "one.sided"), tol=0.00001, ...) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
# The statistics are standardized Brownian motion in [0,1] time interval
# The times are the proportion of information available.
# initialize the correlation matrix
  if (nlook == 1) {
    if (alternative=="one.sided") {
      zz <- qnorm(1-eprob)
    } else {
      zz <- qnorm(1-eprob/2)
    }
  } else {
    corr <- diag(1,nlook)
    for(i in 2:nlook) {
      for(j in 1:(i-1)) {
        corr[i,j] <- min(ifrac[i],ifrac[j])/sqrt(ifrac[i]*ifrac[j])
      }
    }
# The boundary is of the form constant/(ifrac^delta) where
# delta =0 for Pocock and =0.5 for O'Brien-Fleming boundary
# Calculate the constant for a given boundary crossing probability
    bdryden <- ifrac^delta
# iter <- 1
    if (alternative=="one.sided") {
      zlo <- qnorm(1-eprob)
      zhi <- qnorm(1-eprob/nlook)
      alo <- 1 - pmvnorm(upper=zlo/bdryden, corr=corr, ...)
      ahi <- 1 - pmvnorm(upper=zhi/bdryden, corr=corr, ...)
      zz <- zlo + (alo-eprob)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(upper=zz/bdryden, corr=corr, ...)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      while(abs(aa-eprob) > tol) {
# iter <- iter + 1
        if (aa>eprob) {
          zlo <- zz
          alo <- aa
        } else {
          zhi <- zz
          ahi <- aa
        }
        zz <- zlo + (alo-eprob)*(zhi-zlo)/(alo-ahi)
        aa <- 1 - pmvnorm(upper=zz/bdryden, corr=corr, ...)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      }
    } else {
      zlo <- qnorm(1-eprob/2)
      zhi <- qnorm(1-eprob/(2*nlook))
      alo <- 1 - pmvnorm(lower=-zlo/bdryden, upper=zlo/bdryden, corr=corr, ...)
      ahi <- 1 - pmvnorm(lower=-zhi/bdryden, upper=zhi/bdryden, corr=corr, ...)
      zz <- zlo + (alo-eprob)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(lower=-zz/bdryden, upper=zz/bdryden, corr=corr, ...)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      while(abs(aa-eprob) > tol) {
# iter <- iter + 1
        if (aa>eprob) {
          zlo <- zz
          alo <- aa
        } else {
          zhi <- zz
          ahi <- aa
        }
        zz <- zlo + (alo-eprob)*(zhi-zlo)/(alo-ahi)
        aa <- 1 - pmvnorm(lower=-zz/bdryden, upper=zz/bdryden, corr=corr, ...)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      }
    }
  }
  zz
}

###############################################################################

gsd.drift.efficacy <- function(ifrac, delta.eb, sig.level=0.05, pow=0.8, alternative=c("two.sided", "one.sided"), tol=0.00001, ...) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
  corr <- diag(1,nlook)
  for(i in 2:nlook) {
    for(j in 1:(i-1)) {
      corr[i,j] <- min(ifrac[i],ifrac[j])/sqrt(ifrac[i]*ifrac[j])
    }
  }
  ebden <- ifrac^delta.eb
  z0eb <- gsd.bdryconstant(ifrac, sig.level, delta.eb, alternative, tol)
  effbdry <- z0eb/ebden
# Under the alternative the Brownian motion has drift mu.
# So the standardized statistics have means mu*sqrt(ifrac)
# Calculate drift for given power.
  if (alternative=="one.sided") {
    drift0 <- qnorm(pow) + qnorm(1-sig.level)
  } else {
    drift0 <- qnorm(pow) + qnorm(1-sig.level/2)
  }
  driftlo <- drift0*0.8
  powlo <- 1 - pmvnorm(upper=effbdry-driftlo*sqrt(ifrac), corr=corr, ...)
  drifthi <- drift0*1.25
  powhi <- 1 - pmvnorm(upper=effbdry-drifthi*sqrt(ifrac), corr=corr, ...)
  drift0 <- driftlo + (pow-powlo)*(drifthi-driftlo)/(powhi-powlo)
  pow0 <- 1 - pmvnorm(upper=effbdry-drift0*sqrt(ifrac), corr=corr, ...)
#  iter <- 1
#  cat(paste("iteration =",iter,"; drift =",drift0,"; pow =",pow0,"\n"))
#  while (abs(pow0-pow)>tol & iter<50) {
  while (abs(pow0-pow)>tol) {
    if (pow0>pow) {
      drifthi <- drift0
      powhi <- pow0
    } else {
      driftlo <- drift0
      powlo <- pow0
    }
    drift0 <- driftlo + (pow-powlo)*(drifthi-driftlo)/(powhi-powlo)
    pow0 <- 1 - pmvnorm(upper=effbdry-drift0*sqrt(ifrac), corr=corr, ...)
#    iter <- iter+1
#    cat(paste("iteration =",iter,"; drift =",drift0,"; pow =",pow0,"\n"))
  }
  attributes(drift0) <- NULL
  list("effbdry"=effbdry, "drift0"=drift0)
}

###############################################################################

gsd.drift.both <- function(ifrac, delta.eb, delta.fb, sig.level=0.05, pow=0.8, alternative=c("two.sided", "one.sided"), tol=0.00001, ...) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
  corr <- diag(1,nlook)
  for(i in 2:nlook) {
    for(j in 1:(i-1)) {
      corr[i,j] <- min(ifrac[i],ifrac[j])/sqrt(ifrac[i]*ifrac[j])
    }
  }
  ebden <- ifrac^delta.eb
  z0eb <- gsd.bdryconstant(ifrac, sig.level, delta.eb, alternative, tol)
  fbden <- ifrac^delta.fb
  z0fb <- gsd.bdryconstant(ifrac, 1-pow, delta.fb, alternative="one.sided", tol)
# Under the alternative the Brownian motion has drift mu.
# So the standardized statistics have means mu*sqrt(ifrac)
# drift is initiated at z0eb + z0fb
# calculate the power for the initial boundary
  lowerbdry <- -z0fb/fbden
# for two sided test need futility bdry above 0
  if (alternative=="two.sided") {
    lowerbdry[lowerbdry + (z0eb + z0fb)*sqrt(ifrac) < 0] <- -Inf
  }
  upperbdry <- z0eb/ebden - (z0eb + z0fb)*sqrt(ifrac)
  pow0 <- 1 - pnorm(upperbdry[1])
  for(i in 2:nlook) {
    pow0 <- pow0 + pmvnorm(lower=c(lowerbdry[1:(i-1)], upperbdry[i]), upper=c(upperbdry[1:(i-1)], Inf), corr=corr[1:i,1:i], ...)
  }
# adjust boundary constants if pow0 is 1% different from specified power
  sfaclo <- 0.95
  lowerbdry <- -z0fb*sfaclo/fbden
  if (alternative=="two.sided") {
    lowerbdry[lowerbdry + (z0eb + z0fb*sfaclo)*sqrt(ifrac) < 0] <- -Inf
  }
  upperbdry <- z0eb/ebden - (z0eb + z0fb*sfaclo)*sqrt(ifrac)
  powlo <- 1 - pnorm(upperbdry[1])
  for(i in 2:nlook) {
    powlo <- powlo + pmvnorm(lower=c(lowerbdry[1:(i-1)], upperbdry[i]), upper=c(upperbdry[1:(i-1)], Inf), corr=corr[1:i,1:i], ...)
  }
  sfachi <- 1.05
  lowerbdry <- -z0fb*sfachi/fbden
  if (alternative=="two.sided") {
    lowerbdry[lowerbdry + (z0eb + z0fb*sfachi)*sqrt(ifrac) < 0] <- -Inf
  }
  upperbdry <- z0eb/ebden - (z0eb + z0fb*sfachi)*sqrt(ifrac)
  powhi <- 1 - pnorm(upperbdry[1])
  for(i in 2:nlook) {
    powhi <- powhi + pmvnorm(lower=c(lowerbdry[1:(i-1)], upperbdry[i]), upper=c(upperbdry[1:(i-1)], Inf), corr=corr[1:i,1:i], ...)
  }
  sfac0 <- iter <- 1
  while (abs(pow0-pow)/pow > tol & iter<10) {
    iter <- iter+1
    sfac0 <- sfaclo + (sfachi-sfaclo)*(pow - powlo)/(powhi-powlo)
    lowerbdry <- -z0fb*sfac0/fbden
    if (alternative=="two.sided") {
      lowerbdry[lowerbdry + (z0eb + z0fb*sfac0)*sqrt(ifrac) < 0] <- -Inf
    }
    upperbdry <- z0eb/ebden - (z0eb + z0fb*sfac0)*sqrt(ifrac)
    pow0 <- 1 - pnorm(upperbdry[1])
    for(i in 2:nlook) {
      pow0 <- pow0 + pmvnorm(lower=c(lowerbdry[1:(i-1)], upperbdry[i]), upper=c(upperbdry[1:(i-1)], Inf), corr=corr[1:i,1:i], ...)
    }
    if (pow0>pow) {
      sfachi <- sfac0
      powhi <- pow0
    } else {
      sfaclo <- sfac0
      powlo <- pow0
    }
  }
  drift0 <- z0eb + z0fb*sfac0
  effbdry <- z0eb/ebden
  futbdry <- -z0fb*sfac0/fbden + drift0*sqrt(ifrac)
  if (alternative=="two.sided") futbdry[futbdry < 0] <- -Inf
  attributes(drift0) <- NULL
  list("effbdry"=effbdry, "futbdry"=futbdry, "drift0"=drift0)
}

###############################################################################

gsd.drift <- function(ifrac, sig.level = 0.05, pow = 0.8, delta.eb = 0.5, delta.fb = NULL, alternative = c("two.sided", "one.sided"), tol=0.00001, ...) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
  futility <- ifelse(is.null(delta.fb) | missing(delta.fb), FALSE, TRUE)
# check that the last value of information fraction is 1
  if (ifrac[nlook] != 1) stop("last information fraction (ifrac) value should be 1")
  if (nlook == 1) {
    out <- list()
    if (alternative=="one.sided") {
      out$effbdry <- qnorm(1-sig.level)
    } else {
      out$effbdry <- qnorm(1-sig.level/2)
    }
    out$drift0 <- out$effbdry + qnorm(pow)
  } else {
    if (!all(diff(ifrac) > 0)) stop("information fraction (ifrac) values should be increasing")
    if (ifrac[1] <= 0) stop("information fraction (ifrac) values should be positive")
    if (futility) {
      out <- gsd.drift.both(ifrac, delta.eb, delta.fb, sig.level, pow, alternative, tol)
    } else {
      out <- gsd.drift.efficacy(ifrac, delta.eb, sig.level, pow, alternative, tol)
    }
  }
  list("ifrac"=ifrac, "sig.level"=sig.level, "power"=pow, "alternative"=alternative, "delta.eb"=delta.eb, "effbdry"=out$effbdry, "delta.fb"=delta.fb, "futbdry"=out$futbdry, "drift"=out$drift0)
}

# drift parameter theta
# binomial: n (per arm) = theta^2 * 2*pbar*(1-pbar)/(pC -pE)^2  (pooled)
#                  theta^2 * (pC*(1-pC)+pE*(1-pE))/(pC -pE)^2 (unpooled)
#                  CPS correction inflates this number
# normal:   n (per arm) = theta^2 * 2*sigma^2/(muC-muE)^2
# survival: d (total) = theta^2 * 4/(log(haz-ratio))^2
#           Convert number of events d to sample size n

gsdesign.binomial <- function(ifrac, pC, pE, sig.level=0.05, power=0.8,
                              delta.eb = 0.5, delta.fb = NULL, alternative =
                              c("two.sided", "one.sided"), pooled.variance =
                              FALSE, CPS = TRUE, tol = 0.00001, ...) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  if (pooled.variance) {
    pbar <- (pC + pE)/2
    n <- 2 * pbar * (1 - pbar) * (drift.out$drift/(pC - pE))^2
  } else {
    n <- (pC * (1 - pC) + pE * (1 - pE)) * (drift.out$drift/(pC - pE))^2
  }
  if (CPS) {
    A <- n * (pC-pE)^2
    n <- n * {{1+sqrt(1+4*abs(pC-pE)/A)}/2}^2
  }
  out <- drift.out
  out$drift <- NULL
  out$pC <- pC
  out$pE <- pE
  out$outcome <- "binary"
  out$sample.size <- n
  class(out) <- "gsdesign"
  out
}

gsdesign.normal <- function(ifrac, delta, sd=1, sig.level=0.05, power=0.8,
                            delta.eb = 0.5, delta.fb = NULL, alternative =
                            c("two.sided", "one.sided"), tol=0.00001, ...) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  n <- 2*(drift.out$drift*sd/delta)^2
  out <- drift.out
  out$drift <- NULL
  out$delta <- delta
  out$sd <- sd
  out$outcome <- "normal"
  out$sample.size <- n
  class(out) <- "gsdesign"
  out
}

gsdesign.survival <- function(ifrac, haz.ratio, sig.level = 0.05, power = 0.8,
                              delta.eb = 0.5, delta.fb = NULL, alternative = 
                              c("two.sided", "one.sided"), tol=0.00001, ...) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  d <- 4*(drift.out$drift/log(haz.ratio))^2
  out <- drift.out
  out$drift <- NULL
  out$haz.ratio <- haz.ratio
  out$outcome <- "survival"
  out$num.events <- d
  class(out) <- "gsdesign"
  out
}

print.gsdesign <- function(x, ...) {
  if (class(x) != "gsdesign") stop("input shoud be a gsdesign class object")
  cat("\n Group sequential design for comparing", x$outcome, "data with ")
  switch(match(x$outcome, c("binary", "normal", "survival")),
         cat("rates  pC =", x$pC,", pE =", x$pE, "\n"),
         cat("delta =", x$delta, ", sd =", x$sd, "\n"),
         cat("hazard ratio =", x$haz.ratio, "\n"))
  cat("   power family of boundary; 0 (Pocock) to 0.5 (O'Brien-Fleming) \n\n")

  switch(match(x$outcome, c("binary", "normal", "survival")),
         cat("  sample size (per arm) =", x$sample.size, "\n"),
         cat("  sample size (per arm) =", x$sample.size, "\n"),
         cat(" total number of events =", x$num.events, "\n"))
  cat("   information fraction =", format(round(x$ifrac, 3), digits=3), "\n")
  cat("      efficacy boundary =", round(x$effbdry, 3), paste("(power = ", x$delta.eb, ")", sep=""), "\n")
  if (!is.null(x$futbdry)) {
    cat("      futility boundary =", round(x$futbdry, 3), paste("(power = ", x$delta.fb,")",sep=""), "\n")
  }
  cat("              sig.level =", x$sig.level, "\n")
  cat("                  power =", x$power, "\n")
  cat("            alternative =", x$alternative, "\n\n")
  invisible(x)
}
