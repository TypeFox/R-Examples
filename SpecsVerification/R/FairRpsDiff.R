################################
#
# ANALYZE DIFFERENCE IN THE FAIR RANKED PROBABILITY SCORE BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATION
#
# ens     ... ensemble to be tested, ens[i,j] is the number of ensemble members that predict category j at time i
# ens.ref ... reference forecast ensemble, same dimension as ens
# obs     ... observations matrix, same dimension as ens, obs[i,j] = 1 if category j occured at time i, 0 otherwise
#
# return value: a list with elements
#     * rps.diff ... the fair Ranked Probability Score difference
#     * sampling.quantiles ... if `probs` were defined, the corresponding
#                              quantiles of the sampling distribution 
#                              of `rps.diff`
#
################################
FairRpsDiff <- function(ens, ens.ref, obs, probs=NA) {

  stopifnot(all(dim(ens) == dim(obs)))
  stopifnot(all(dim(ens.ref) == dim(obs)))


  # calculate fair RPS score differences
  rps.ens <- FairRps(ens, obs)
  rps.ref <- FairRps(ens.ref, obs)
  rps.diff <- rps.ref - rps.ens
  N <- nrow(obs) - sum(is.na(rps.diff))
  mean.rps.diff <- mean(rps.diff, na.rm=TRUE)

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs))) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd(rps.diff) / sqrt(N) + mean.rps.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  if (N > 1) {
    p.value <- 1-pt(mean.rps.diff / sd(rps.diff) * sqrt(N), df=N-1)
  } else {
    p.value <- NA
  }

  #return
  list(rps.diff=mean.rps.diff, sampling.quantiles=cis, p.value=p.value)
}

