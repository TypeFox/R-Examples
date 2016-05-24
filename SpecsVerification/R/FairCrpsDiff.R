################################
#
# ANALYZE DIFFERENCE IN THE FAIR CRPS BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATIONS
#
# ens     ... the ensemble (matrix of dimension N*K)
# ens.ref ... the reference ensemble (matrix of dimension N*K.ref)
# obs     ... observations (vector of length N)
# probs   ... quantiles of the sampling distribution
#
################################
FairCrpsDiff <- function(ens, ens.ref, obs, probs=NA) {


  # pre-processing
  l <- Preprocess(ens=ens, ens.ref=ens.ref, obs=obs)
  ens <- l[["ens"]]
  ens.ref <- l[["ens.ref"]]
  obs <- l[["obs"]]

  N <- length(obs)

  # calculate fair crps difference
  crps.ens <- FairCrps(ens, obs)
  crps.ref <- FairCrps(ens.ref, obs)
  crps.diff <- crps.ref - crps.ens
  mean.crps.diff <- mean(crps.diff, na.rm=TRUE)
  sd.crps.diff <- sd(crps.diff, na.rm=TRUE)

  # update N, accounting for NA's
  N <- N - sum(is.na(crps.diff))

  # if all scores are NA, return NA
  if (all(is.na(crps.diff))) {
    return(list(crps.diff=NA, sampling.quantiles=probs*NA, p.value=NA))
  }

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs)) && N > 1) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd.crps.diff / sqrt(N) + mean.crps.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  p.value <- ifelse(N>1, 
    1-pt(mean.crps.diff / sd.crps.diff * sqrt(N), df=N-1), 
    NA)

  #return
  list(crps.diff=mean.crps.diff, sampling.quantiles=cis, p.value=p.value)
}

