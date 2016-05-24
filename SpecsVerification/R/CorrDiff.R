################################
#
# ANALYZE DIFFERENCE IN THE CORRELATION COEFFICIENT BETWEEN TWO ENSEMBLE MEANS
# FOR THE SAME OBSERVATIONS
#
# ens     ... the ensemble (matrix of dimension N*K)
# ens.ref ... the reference ensemble (matrix of dimension N*K.ref)
# obs     ... observations (vector of length N)
# sign.level ... significance level of the confidence interval
#
################################
CorrDiff <- function(ens, ens.ref, obs, sign.level=0.05) {

  # preprocess
  l <- Preprocess(ens=ens, ens.ref=ens.ref, obs=obs)
  ens <- l[["ens"]]
  ens.ref <- l[["ens.ref"]]
  obs <- l[["obs"]]

  N <- length(obs)

  # calculate correlation coefficients and their confidence intervals
  cc.ens <- Corr(ens, obs, probs=c(0.5*sign.level, 1-0.5*sign.level))[1:3]
  cc.ref <- Corr(ens.ref, obs, probs=c(0.5*sign.level, 1-0.5*sign.level))[1:3]

  # calculate correlation difference
  cc.diff <- cc.ens[1] - cc.ref[1]

  # auxiliary quantities
  r12 <- cc.ens[1]
  r13 <- cc.ref[1]
  r23 <- cor(rowMeans(ens, na.rm=TRUE), rowMeans(ens.ref, na.rm=TRUE), 
             use="pairwise.complete.obs")

  # confidence interval, according to zou 2007, example 2
  if (sign.level <= 0 | sign.level >= 1) {
    sign.level <- NA
  }
  if (is.na(sign.level)) {
    L <- U <- NA
  } else {
    # individual confidence limits of cc.ens and cc.ref
    l1 <- cc.ens[2]
    u1 <- cc.ens[3]
    l2 <- cc.ref[2]
    u2 <- cc.ref[3]
  
    # correlation between the two corcoefs r12 and r13
    c.12.13 <- ((r23 - 0.5 * r12 * r13) * (1 - r12*r12 - 
      r13*r13 - r23*r23) + r23*r23*r23) / ((1 - r12*r12) * 
      (1 - r13*r13))
    # lower confidence limit
    L <- r12 - r13 - sqrt((r12-l1)^2 + (u2 - r13)^2 - 
      2*c.12.13*(r12-l1)*(u2-r13))
    # upper confidence limit
    U <- r12 - r13 + sqrt((u1 - r12)^2 + (r13-l2)^2 - 
      2*c.12.13*(u1-r12)*(r13-l2))
  }

  # N minus number of NA rows
  N <- length(obs) - sum(is.na(rowSums(ens, na.rm=TRUE) + 
       rowSums(ens.ref, na.rm=TRUE) + obs))

  # p value of one-sided test for equality of dependent correlation
  # coefficients (steiger 1980 Eq 7)
  if (N < 4) {
    p.value <- NA
  } else {
    R <- (1-r12*r12-r13*r13-r23*r23) + 2*r12*r13*r23
    t <- (r12 - r13) * sqrt((N-1)*(1+r23) / (2 * ((N-1)/(N-3)) * 
      R + 0.25*(r12+r13)^2 * (1-r23)^3))
    p.value <- 1 - pt(t, df=N-3)
  }


  #return
  ret <- c(cc.diff, L, U, p.value)
  names(ret) <- c("corr.diff", c("L","U"), "p.value")
  return(ret)
}


