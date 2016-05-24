################################
#
# ANALYZE DIFFERENCE IN THE CRPS BETWEEN TWO DRESSED ENSEMBLES
# FOR THE SAME VERIFICATION
#
# ens     ... dressed ensemble (object of class `dressed.ens`)
# ens.ref ... dressed reference ensemble (object of class `dressed.ens`)
# obs     ... verifications (vector of length N)
# probs   ... quantiles of the sampling distribution
#
################################
DressCrpsDiff <- function(dressed.ens, dressed.ens.ref, obs, probs=NA) {

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
  crps.ens <- DressCrps(dressed.ens, obs)
  crps.ref <- DressCrps(dressed.ens.ref, obs)
  crps.diff <- crps.ref - crps.ens
  mean.crps.diff <- mean(crps.diff)

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs))) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd(crps.diff) / sqrt(N) + mean.crps.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  p.value <- 1-pt(mean.crps.diff / sd(crps.diff) * sqrt(N), df=N-1)

  #return
  list(crps.diff=mean.crps.diff, sampling.quantiles=cis, p.value=p.value)
}

