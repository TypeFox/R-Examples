################################
#
# ANALYZE DIFFERENCE IN THE IGNORANCE BETWEEN TWO DRESSED ENSEMBLES
# FOR THE SAME VERIFICATION
#
# ens     ... dressed ensemble (object of class `dressed.ens`)
# ens.ref ... dressed reference ensemble (object of class `dressed.ens`)
# obs     ... verifications (vector of length N)
# probs   ... quantiles of the sampling distribution
#
################################
DressIgnDiff <- function(dressed.ens, dressed.ens.ref, obs, probs=NA) {

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
  ign.diff <- ign.ref - ign.ens
  mean.ign.diff <- mean(ign.diff)

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs))) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd(ign.diff) / sqrt(N) + mean.ign.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  p.value <- 1-pt(mean.ign.diff / sd(ign.diff) * sqrt(N), df=N-1)

  #return
  list(ign.diff=mean.ign.diff, sampling.quantiles=cis, p.value=p.value)
}

