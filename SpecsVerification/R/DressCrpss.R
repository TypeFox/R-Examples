################################
#
# ANALYZE CRPS IMPROVEMENT BETWEEN TWO DRESSED ENSEMBLES FOR THE SAME
# VERIFICATION BY THE SKILL SCORE 1 - S / S.ref
#
# ens     ... dressed ensemble (object of class `dressed.ens`)
# ens.ref ... dressed reference ensemble (object of class `dressed.ens`)
# obs     ... verifications (vector of length N)
#
################################
DressCrpss <- function(dressed.ens, dressed.ens.ref, obs) {

  # sanity checks
  if (class(obs) == "data.frame") {
    obs <- c(as.matrix(obs))
  }
  stopifnot(class(dressed.ens)=="dressed.ens", 
            class(dressed.ens.ref)=="dressed.ens")
  stopifnot(is.vector(obs), length(obs) > 1)

  ens <- dressed.ens[["ens"]]
  ens.ref <- dressed.ens.ref[["ens"]]

  stopifnot(nrow(ens)==length(obs), nrow(ens.ref) == length(obs))

  N <- length(obs)
  K <- ncol(ens)
  K.ref <- ncol(ens.ref)
  obs <- matrix(obs, ncol=1)

  K <- ncol(ens)
  K.ref <- ncol(ens.ref)

  # calculate crpss
  crps.ens <- DressCrps(dressed.ens, obs)
  crps.ref <- DressCrps(dressed.ens.ref, obs)
  crpss <- 1 - mean(crps.ens) / mean(crps.ref)
  crpss.sigma <- 1 / sqrt(N) * sqrt( var(crps.ens) / mean(crps.ref)^2 + 
         var(crps.ref) * mean(crps.ens)^2 / mean(crps.ref)^4 - 
         2 * cov(crps.ens, crps.ref) * mean(crps.ens) / mean(crps.ref)^3)

  #return
  list(crpss=crpss, crpss.sigma=crpss.sigma)
}

