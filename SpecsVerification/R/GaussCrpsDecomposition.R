################################
#
# CRPS DECOMPOSITION OF GAUSSIAN FORECASTS WITH MEANS `mean`, 
# STANDARD DEVIATIONS `sd`, FOR VERIFYING OBSERVATIONS `obs`
#
# mean ... vector of forecast means (or single value if mean is constant)
# sd   ... vector of forecast standard deviations (or single value if sd is
#          constant)
# obs  ... vector of observations
# probs  ... which quantiles of the sampling distributions should be
#            estimated by bootstrapping?
# n.boot ... number of bootstrap replicates 
#
################################

GaussCrpsDecomposition <- function(mean,sd,obs,probs=NA,n.boot=0) {

  # CRPS divergence function between N(m,s^2) and N(n,t^2)
  # if s=t=0, then |m-n| is returned
  crps.div <- function(m, s, n, t) {
    sq.st <- sqrt(s * s + t * t)
    z <- (m - n) / sq.st
    div <- sq.st * (z * (2 * pnorm(z) - 1) + 
           2 * dnorm(z)) - (s + t) / sqrt(pi)
    i.zero <- which(sq.st==0)
    div[i.zero] <- (abs(m-n))[i.zero]
    return(div)
  }

  # transform data frames and matrices to vectors
  mean <- c(as.matrix(mean))
  sd <- c(as.matrix(sd))
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

  # check if all vectors are of equal length
  len <- c(length(mean), length(sd), length(obs))
  if (length(unique(len)) > 1) {
    stop("input vectors have different lengths")
  }

  # replace NA, NaN, Inf by NA
  mean[ !is.finite(mean) ] <- NA
  sd[ !is.finite(sd) ] <- NA
  obs[ !is.finite(obs) ] <- NA

  # calibration by minimum CRPS estimation
  crps.opt <- function(par) {
    with(as.list(par), {
      m.cur <- a + b * mean
      s.cur <- c + d * d * sd
      mean(crps.div(m.cur, s.cur, obs, 0), na.rm=TRUE)
    })
  }
  # start from no calibration and optimize using BFGS
  par <- c(a=0, b=1, c=0, d=1) 
  par <- optim(par=par, fn=crps.opt, method="BFGS")[["par"]]
  mean.cal <- par[1] + par[2] * mean
  sd.cal <- par[3] + par[4]^2 * sd

  # climatological parameters
  mean.clim <- rep(mean(obs, na.rm=TRUE), N)
  sd.clim <- rep(sd(obs, na.rm=TRUE), N)

  # the decomposition
  rel <- mean(crps.div(mean, sd, mean.cal, sd.cal), na.rm=TRUE)
  res <- mean(crps.div(mean.cal, sd.cal, mean.clim, sd.clim), na.rm=TRUE)
  unc <- sd.clim[1] / sqrt(pi)
  crps <- mean(crps.div(mean, sd, obs, 0), na.rm=TRUE)
  crps.clim <- mean(crps.div(mean.clim, sd.clim, obs, 0), na.rm=TRUE)
  crps.cal <- mean(crps.div(mean.cal, sd.cal, obs, 0), na.rm=TRUE)


  # quantiles of the sampling distribution 
  boot.qntls <- NA
  if (!any(is.na(probs)) & n.boot > 1) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    boot <- t(replicate(n.boot, {

      inds <- sample(1:N, N, replace=TRUE)
      bmean <- mean[inds]
      bsd <- sd[inds]
      bobs <- obs[inds]
    
      # calibration by minimum CRPS estimation
      crps.opt <- function(par) {
        with(as.list(par), {
          m.cur <- a + b * bmean
          s.cur <- c + d * d * bsd
          mean(crps.div(m.cur, s.cur, bobs, 0), na.rm=TRUE)
        })
      }
      # start from no calibration and optimize using BFGS
      par <- c(a=0, b=1, c=0, d=1) 
      par <- optim(par=par, fn=crps.opt, method="BFGS")[["par"]]
      mean.cal <- par[1] + par[2] * bmean
      sd.cal <- par[3] + par[4]^2 * bsd
    
      # climatological parameters
      mean.clim <- rep(mean(bobs, na.rm=TRUE), N)
      sd.clim <- rep(sd(bobs, na.rm=TRUE), N)
    
      # the decomposition
      rel <- mean(crps.div(bmean, bsd, mean.cal, sd.cal), na.rm=TRUE)
      res <- mean(crps.div(mean.cal, sd.cal, mean.clim, sd.clim), na.rm=TRUE)
      unc <- sd.clim[1] / sqrt(pi)
      crps <- mean(crps.div(bmean, bsd, bobs, 0), na.rm=TRUE)
      crps.clim <- mean(crps.div(mean.clim, sd.clim, bobs, 0), na.rm=TRUE)
      crps.cal <- mean(crps.div(mean.cal, sd.cal, bobs, 0), na.rm=TRUE)
      c(REL=rel, RES=res, UNC=unc, CRPS=crps, CRPS.clim=crps.clim, CRPS.cal=crps.cal)
    }))# end of bootstrap

    boot.qntls <- apply(boot, 2, quantile, probs=probs)
    rownames(boot.qntls) <- paste(probs)
  }


  # return
  list(REL=rel, RES=res, UNC=unc, CRPS=crps, ngr.pars=par, 
       CRPS.clim=crps.clim, CRPS.cal=crps.cal, BOOT.quantiles=boot.qntls)
}

