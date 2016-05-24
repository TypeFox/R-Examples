################################
#
# ANALYZE DIFFERENCE IN THE CRPS BETWEEN TWO FORECASTS ISSUED AS NORMAL
# DISTRIBUTIONS FOR THE SAME OBSERVATIONS
#
# mean     ... forecast mean (vector of length N)
# sd       ... forecast standard deviation (vector of length N)
# mean.ref ... mean of the reference forecast (vector of length N)
# sd.ref   ... standard deviation of the reference forecast (vector of length N)
# obs      ... observations (vector of length N)
# probs    ... quantiles of the sampling distribution
#
################################
GaussCrpsDiff <- function(mean, sd, mean.ref, sd.ref, obs, probs=NA) {

  # transform data frames and matrices to vectors
  mean <- c(as.matrix(mean))
  mean.ref <- c(as.matrix(mean.ref))
  sd <- c(as.matrix(sd))
  sd.ref <- c(as.matrix(sd.ref))
  obs <- c(as.matrix(obs))

  # 
  N <- length(obs)

  # if any forecast vector is of length one, expand them to length N
  if (length(mean) == 1) {
    mean <- rep(mean, N)
  }
  if (length(sd) == 1) {
    sd <- rep(sd, N)
  }
  if (length(mean.ref) == 1) {
    mean.ref <- rep(mean.ref, N)
  }
  if (length(sd.ref) == 1) {
    sd.ref <- rep(sd.ref, N)
  }
  

  # check if all vectors are of equal length
  len <- c(length(mean), length(sd), length(obs), length(mean.ref), length(sd.ref))
  if (length(unique(len)) > 1) {
    stop("vectors have different lengths")
  }

  # calculate crps difference
  crps <- GaussCrps(mean, sd, obs)
  crps.ref <- GaussCrps(mean.ref, sd.ref, obs)
  crps.diff <- crps.ref - crps
  mean.crps.diff <- mean(crps.diff, na.rm=TRUE)

  # reduce N if there are NA's
  N <- N - sum(is.na(crps.diff))

  # quantiles of the sampling distribution 
  cis <- NA
  if (!any(is.na(probs))) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    cis <- qt(probs, df=N-1) * sd(crps.diff, na.rm=TRUE) / sqrt(N) + mean.crps.diff
    names(cis) <- paste(probs)
  }

  # p value of paired one-sided t test for positive score difference
  p.value <- 1-pt(mean.crps.diff / sd(crps.diff) * sqrt(N), df=N-1)

  #return
  list(crps.diff=mean.crps.diff, sampling.quantiles=cis, p.value=p.value)
}

