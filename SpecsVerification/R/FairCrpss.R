################################
#
# ANALYZE DIFFERENCE IN THE FAIR CRPS BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATIONS BY SKILL SCORE
#
# ens     ... the ensemble (matrix of dimension N*K)
# ens.ref ... the reference ensemble (matrix of dimension N*K.ref)
# obs     ... observations (vector of length N)
#
################################
FairCrpss <- function(ens, ens.ref, obs) {

  # pre-process
  l <- Preprocess(ens=ens, ens.ref=ens.ref, obs=obs) 
  ens <- l[["ens"]]
  ens.ref <- l[["ens.ref"]]
  obs <- l[["obs"]]

  # calculate individual scores
  crps.ens <- FairCrps(ens, obs)
  crps.ref <- FairCrps(ens.ref, obs)

  # only use pairwise complete scores
  i.na <- !(is.na(crps.ens + crps.ref))
  crps.ens <- crps.ens[i.na]
  crps.ref <- crps.ref[i.na]

  if (!any(i.na)) {
    return(list(bss=NA, bss.sigma=NA))
  }


  # calculate auxiliary quantities
  N <- length(obs)
  m.crps.ens <- mean(crps.ens)
  m.crps.ref <- mean(crps.ref)
  v.crps.ens <- var(crps.ens)
  v.crps.ref <- var(crps.ref)
  cov.crps   <- cov(crps.ens, crps.ref)

  # calculate skill score
  crpss <- 1 - m.crps.ens / m.crps.ref

  # update N
  N <- N - sum(is.na(crps.ens+crps.ref))


  # calculate error propagation standard deviation
  crpss.sigma <- ifelse(N > 1,
         1 / sqrt(N) * sqrt( v.crps.ens / m.crps.ref^2 + 
         v.crps.ref * m.crps.ens^2 / m.crps.ref^4 - 
         2 * cov.crps * m.crps.ens / m.crps.ref^3),
         NA)

  #return
  list(crpss=crpss, crpss.sigma=crpss.sigma)
}

