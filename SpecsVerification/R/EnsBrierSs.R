################################
#
# ANALYZE DIFFERENCE IN THE BRIER SCORE BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME OBSERVATION BY SKILL SCORE 1 - S / S.ref
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
#     * br.ss ... the Brier Skill Score 
#     * br.sigma ... approximate standard deviation of the skill score
#
################################
EnsBrierSs <- function(ens, ens.ref, obs, tau=0.5) {

  # pre-process
  l <- Preprocess(ens=ens, ens.ref=ens.ref, obs=obs) 
  ens <- l[["ens"]]
  ens.ref <- l[["ens.ref"]]
  obs <- l[["obs"]]

  # sanity checks
  stopifnot(length(tau) == 1 | length(tau) == length(obs))


  # calculate Brier skill score 
  br.ens <- EnsBrier(ens, obs, tau)
  br.ref <- EnsBrier(ens.ref, obs, tau)

  # only use pairwise complete scores
  i.na <- !(is.na(br.ens + br.ref))
  br.ens <- br.ens[i.na]
  br.ref <- br.ref[i.na]

  if (!any(i.na)) {
    return(list(bss=NA, bss.sigma=NA))
  }

  # calculate auxiliary quantities
  N <- length(obs)
  m.br.ens <- mean(br.ens)
  m.br.ref <- mean(br.ref)
  v.br.ens <- var(br.ens)
  v.br.ref <- var(br.ref)
  cov.br   <- cov(br.ens, br.ref)

  # calculate skill score
  bss <- 1 - m.br.ens / m.br.ref

  # update N
  N <- N - sum(is.na(br.ens+br.ref))

  # calculate error propagation standard deviation
  bss.sigma <- ifelse(N > 1,
         1 / sqrt(N) * sqrt( v.br.ens / m.br.ref^2 + 
         v.br.ref * m.br.ens^2 / m.br.ref^4 - 
         2 * cov.br * m.br.ens / m.br.ref^3),
         NA)

  #return
  list(bss=bss, bss.sigma=bss.sigma)
}

