################################
#
# ANALYZE IGNORANCE IMPROVEMENT BETWEEN TWO DRESSED ENSEMBLES FOR THE SAME
# VERIFICATION BY THE SKILL SCORE 1 - S / S.ref
#
# ens     ... dressed ensemble (object of class `dressed.ens`)
# ens.ref ... dressed reference ensemble (object of class `dressed.ens`)
# obs     ... verifications (vector of length N)
#
################################
DressIgnSs <- function(dressed.ens, dressed.ens.ref, obs) {

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

  # calculate crps difference
  ign.ens <- DressIgn(dressed.ens, obs)
  ign.ref <- DressIgn(dressed.ens.ref, obs)


  ignss <- 1 - mean(ign.ens) / mean(ign.ref)
  ignss.sigma <- 1 / sqrt(N) * sqrt( var(ign.ens) / mean(ign.ref)^2 + 
         var(ign.ref) * mean(ign.ens)^2 / mean(ign.ref)^4 - 
         2 * cov(ign.ens, ign.ref) * mean(ign.ens) / mean(ign.ref)^3)

  #return
  list(ignss=ignss, ignss.sigma=ignss.sigma)

}

