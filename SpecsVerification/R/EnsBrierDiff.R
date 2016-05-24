################################
#
# ANALYZE DIFFERENCE IN THE BRIER SCORE BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATION
#
# ens     ... ensemble to be tested (matrix of dimension N*K)
# ens.ref ... reference forecast ensemble (matrix of dimension N*K.ref)
# obs     ... observations (vector of length N)
# tau     ... threshold, whose exceedance defines the "event" (scalar, or
#             vector of length N)
#             the default is 0.5, such that ensemble members can be given as
#             event indicators, i.e. 0 or 1
#
# return value: a list with elements
#     * br.diff ... the Brier Score difference
#     * sampling.quantiles ... if `probs` were defined, the corresponding
#                              quantiles of the sampling distribution 
#                              of `br.diff`
#
################################
EnsBrierDiff <- function(ens, ens.ref, obs, tau=0.5, probs=NA) {

  # pre-processing
  l <- Preprocess(ens=ens, ens.ref=ens.ref, obs=obs)
  ens <- l[["ens"]]
  ens.ref <- l[["ens.ref"]]
  obs <- l[["obs"]]

  # sanity checks
  stopifnot(is.numeric(c(ens, ens.ref, obs, tau)))
  stopifnot(length(tau) == 1 | length(tau) == length(obs))

  N <- length(obs)

  # calculate Brier score differences
  br.ens <- EnsBrier(ens, obs, tau)
  br.ref <- EnsBrier(ens.ref, obs, tau)
  br.diff <- br.ref - br.ens

  if (all(is.na(br.diff))) {
    return(list(br.diff=NA, sampling.quantiles=probs*NA, p.value=NA))
  }

  # calculate mean and sd
  mean.br.diff <- mean(br.diff, na.rm=TRUE)
  sd.br.diff <- sd(br.diff, na.rm=TRUE)
  # update N, accounting for NA's
  N <- N - sum(is.na(br.diff))

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs)) && N > 1) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd.br.diff / sqrt(N) + mean.br.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  p.value <- ifelse(N>1, 
    1-pt(mean.br.diff / sd.br.diff * sqrt(N), df=N-1), 
    NA)

  #return
  list(br.diff=mean.br.diff, sampling.quantiles=cis, p.value=p.value)
}

