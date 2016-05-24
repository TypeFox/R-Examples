################################
#
# ANALYZE DIFFERENCE IN THE RANKED PROBABILITY SCORE BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATION BY SKILL SCORE 1 - S / S.ref
#
# ens     ... ensemble to be tested, ens[i,j] is the number of ensemble members that predict category j at time i
# ens.ref ... reference forecast ensemble, same dimension as ens
# obs     ... observations matrix, same dimension as ens, obs[i,j] = 1 if category j occured at time i, 0 otherwise
#
# return value: a list with elements
#     * rpss ... the Ranked Probability Skill Score 
#     * sigma.rpss ... estimated standard deviation of Ranked Probability Skill Score 
#
################################
EnsRpss <- function(ens, ens.ref, obs) {

  stopifnot(all(dim(ens) == dim(obs)))
  stopifnot(all(dim(ens.ref) == dim(obs)))


  # calculate RPS score differences
  rps.ens <- EnsRps(ens, obs)
  rps.ref <- EnsRps(ens.ref, obs)

  # only use cases where both scores are not NA
  i.nna <- which(!is.na(rps.ens + rps.ref))
  rps.ens <- rps.ens[i.nna]
  rps.ref <- rps.ref[i.nna]
  N <- length(i.nna)

  # calculate skill score
  rpss <- 1 - mean(rps.ens) / mean(rps.ref)
  # estimate standard deviation by propagation of uncertainty 
  rpss.sigma <- 1 / sqrt(N) * sqrt( var(rps.ens) / mean(rps.ref)^2 + 
         var(rps.ref) * mean(rps.ens)^2 / mean(rps.ref)^4 - 
         2 * cov(rps.ens, rps.ref) * mean(rps.ens) / mean(rps.ref)^3)

  #return
  list(rpss=rpss, rpss.sigma=rpss.sigma)

}

