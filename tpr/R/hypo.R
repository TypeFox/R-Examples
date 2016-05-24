############################################################################
## Nonparametric test of significance: H0: beta(t) = c
############################################################################
sig.test.int.1.ff <- function(fit, chypo, idx, weight=TRUE, ncut=2) {
  ## fit:    the returned object from tpr
  ## chypo:  the hopothesized value c: H_0: \alpha[idx] == c
  ## idx:    the index of the coefficient to be tested
  ## weight: whether or not use weight = inverse standard error
  ## ncut:   the number of cuts of the interval of interest
  tis <- fit$tis
  nt <- length(tis)
  alpha.j <- fit$alpha[,idx]
  valpha.j <- fit$valpha[,idx]
  if (weight) w <- 1/sqrt(valpha.j)
  else w <- rep(1, length(sqrt(valpha.j)))

  ## getting the influence functions
  nn <- ncol(fit$inflTheta)
  infls <- lapply(fit$inflAlpha, function (x) x[idx,])
  infls.alpha.j <- matrix(unlist(infls), ncol=nt)
  infls.eta <- chypo

  stat <- rep(NA, ncut + 1)
  for (i in 1:(ncut+1)) {
    wi <- w * (tis >= quantile(tis, (i-1)/(ncut + 1)) & tis < quantile(tis, i / (ncut + 1)))
    infls.i <- apply((infls.alpha.j - infls.eta) * t(matrix(wi, nrow=nt, ncol=nn)), 1, function(x) sum(x[-nt] * diff(tis)))
    Delta.i <- (alpha.j - chypo)
    Delta.i <- sum(Delta.i[-nt] * wi[-nt] * diff(tis))  #### observed statistic
    vDelta.i <- var(infls.i) / nn
    stat[i] <- Delta.i / sqrt(vDelta.i)
  }
  stat
}

sig.test.int.ff <- function(fit, chypo=0, idx, weight=TRUE, ncut=2) {
  stat <- NULL
  nidx <- length(idx)
  ncol <- (ncut + 1) * 2
  chypo <- rep(chypo, length=nidx)
  for (i in 1:nidx) {
    ##cat ("current coefficient", idx[i], "\n")
    foo <- sig.test.int.1.ff(fit, chypo[i], idx[i], weight, ncut)
    stat <- rbind(stat, foo)
  }
  pvalue <- 2 * (1 - pnorm(abs(stat)))

  ret <- matrix(NA, nidx, ncol)
  ret[, 2 * (1:(ncut + 1)) - 1] <- stat
  ret[, 2 * (1:(ncut + 1))] <- pvalue
  row.names(ret) <- idx
  colnames(ret) <- paste(c("stat", "pvalue"), rep(1:(ncut + 1), rep(2, (ncut + 1))), sep=".")
  ret
}

##########################################################################
## Goodness-of-fit test for constant fit: Hc: beta(t) = beta
##########################################################################
gof.test.int.1.ff <- function(fit, cfit, idx, weight=TRUE, ncut=1) {
  ## fit:    the returned object from tpr
  ## cfit:   the constant fit object from cst.fit with influence functions available
  ## idx:    the index of the coefficient to be tested
  ## weight: whether or not use weight = inverse standard error
  ## ncut:   the number of cuts of the interval of interest
  tis <- fit$tis
  nt <- length(tis)
  alpha.j <- fit$alpha[,idx]
  valpha.j <- fit$valpha[,idx]
  if (weight) w <- 1/sqrt(valpha.j)
  else w <- rep(1, length(sqrt(valpha.j)))

  eta <- cfit$fit$eta
  
  ## getting the influence functions
  nn <- ncol(fit$inflTheta)
  infls <- lapply(fit$inflAlpha, function (x) x[idx,])
  infls.alpha.j <- matrix(unlist(infls), ncol=nt)
  infls.eta <- cfit$infls

  stat <- rep(NA, ncut + 1)
  for (i in 1:(ncut+1)) {
    wi <- w * (tis >= quantile(tis, (i-1)/(ncut + 1)) & tis < quantile(tis, i / (ncut + 1)))
    infls.i <- apply((infls.alpha.j - infls.eta) * t(matrix(wi, nrow=nt, ncol=nn)), 1, function(x) sum(x[-nt] * diff(tis)))
    Delta.i <- (alpha.j - eta)
    Delta.i <- sum(Delta.i[-nt] * wi[-nt] * diff(tis))  #### observed statistic
    vDelta.i <- var(infls.i) / nn
    stat[i] <- Delta.i / sqrt(vDelta.i)
  }
  stat
}

gof.test.int.ff <- function(fit, cfitList=NULL, idx, weight=TRUE, ncut=2) {
  ## fit:    the returned object from tpr
  ## cfit:   a list of cfit
  ## idx:    the index of the coefficient to be tested
  ## weight: whether or not use weight = inverse standard error
  ## ncut:   the number of cuts of the interval of interest
  stat <- NULL
  nidx <- length(idx)
  ncol <- (ncut + 1) * 2

  if (is.null(cfitList)) cfitList <- cst.fit.ff(fit, idx)
  for (i in 1:nidx) {
    ##cat ("current coefficient", idx[i], "\n")
    foo <- gof.test.int.1.ff(fit, cfitList[[i]], idx[i], weight, ncut)
    stat <- rbind(stat, foo)
  }
  pvalue <- 2 * (1 - pnorm(abs(stat)))

  ret <- matrix(NA, nidx, ncol)
  ret[, 2 * (1:(ncut + 1)) - 1] <- stat
  ret[, 2 * (1:(ncut + 1))] <- pvalue
  row.names(ret) <- idx
  colnames(ret) <- paste(c("stat", "pvalue"), rep(1:(ncut + 1), rep(2, (ncut + 1))), sep=".")
  ret
}



#########################################################################
## Significance or goodness-of-fit test based on bootstrapping
#########################################################################
test.boots.1 <- function(fit, chypo=0, idx, nsim, cfit=NULL) {
  tis <- fit$tis
  nt <- length(tis)
  alpha.j <- fit$alpha[,idx]
  valpha.j <- fit$valpha[,idx]
  w <- 1/sqrt(valpha.j)

  ## getting the influence functions
  nn <- ncol(fit$inflTheta)
  infls <- lapply(fit$inflAlpha, function (x) x[idx,])
  infls.alpha.j <- matrix(unlist(infls), ncol=nt)

  eta <- ifelse(is.null(cfit), chypo, cfit$fit$eta)
  infls.eta <- ifelse(is.null(cfit), chypo, cfit$infls)

  infl <- infls.alpha.j - infls.eta
  sig <- apply(infl, 2, sd, na.rm = TRUE) / sqrt(nn)
  
  obs <- (alpha.j - eta) / sig

  sim <- .Call("bootsSample_rap", infl, sig, as.integer(nsim), PACKAGE="tpr")

  supobs <- max(abs(obs))
  supsim <- apply(abs(sim), 2, max)
  suppval <- sum(supsim >= supobs) / nsim

  absobs <- sum(abs(obs))
  abssim <- apply(abs(sim), 2, sum)
  abspval <- sum(abssim >= absobs) / nsim

  s2obs <- sum(obs^2)
  s2sim <- apply(sim^2, 2, sum)
  s2pval <- sum(s2sim >= s2obs) / nsim

  intobs <- sum(obs)
  intsim <- apply(sim, 2, sum)
  intpval <- sum(abs(intsim) >= abs(intobs)) / nsim
  
  val <- list(tis = tis, obs = obs, sim = sim,
              supobs = supobs, supsim = supsim, suppval = suppval,
              absobs = absobs, abssim = abssim, abspval = abspval,
              s2obs = s2obs, s2sim = s2sim, s2pval = s2pval,
              intobs = intobs, intsim = intsim, intpval = intpval)
  #print(c(suppval, abspval, s2pval, intpval))
  val
}


sig.test.boots.ff <- function(fit, chypo=0, idx, nsim=1000, plot=FALSE) {
  val <- NULL
  for (i in idx) {
    bar <- test.boots.1(fit=fit, chypo=chypo, idx=i, nsim=nsim)
    val <- rbind(val, c(bar$suppval, bar$abspval, bar$s2pval, bar$intpval,
                        bar$supobs, bar$absobs, bar$s2obs, bar$intobs))
    if (plot) plotGof(bar, 25)
  }
  colnames(val) <- c("sup.pvalue", "abs.pvalue", "square.pvalue", "int.pvalue",
                     "sup.obs", "abs.obs", "square.obs", "int.obs")
  row.names(val) <- idx
  val
}


gof.test.boots.ff <- function(fit, cfitList=NULL, idx, nsim=1000, plot=FALSE) {
  val <- NULL
  if (is.null(cfitList)) cfitList <- cst.fit.ff(fit, idx)
  for (i in 1:length(idx)) {
    bar <- test.boots.1(fit=fit, idx=idx[i], nsim=nsim, cfit=cfitList[[i]])
    val <- rbind(val, c(bar$suppval, bar$abspval, bar$s2pval, bar$intpval,
                        bar$supobs, bar$absobs, bar$s2obs, bar$intobs))
    
    if (plot) plotGof(bar, 25)
  }
  colnames(val) <- c("sup.pvalue", "abs.pvalue", "square.pvalue", "int.pvalue",
                     "sup.obs", "abs.obs", "square.obs", "int.obs")
  row.names(val) <- idx
  val
}


#######################################################################
## goodness-of-fit test for partly functional models
#######################################################################
gof.test.boots.pf <- gofTest <- function(fit1, fit2, nsim, p=NULL, q=1) {
#### fit1 is the H0 model, reduced (constant coef) model
#### fit2 is the full model
#### nsim is the number of bootstrap sample
#### p is the position of the time-varying estimation in fit2
#### q is the position of the time-independent estimator in fit1
  infl1 <- fit1$inflTheta[q,]
  n <- length(infl1)
  if (is.null(p))  p <- ncol(fit2$alpha)
  infl2 <- sapply(fit2$inflAlpha, function(x) x[p,])
  infl <- infl2 - infl1
  sig <- apply(infl, 2, sd, na.rm = TRUE) / sqrt(n)
  obs <- (fit2$alpha[,p] - fit1$theta[q]) / sig
  sim <- .Call("bootsSample_rap", infl, sig, as.integer(nsim), PACKAGE="tpr")
  supobs <- max(abs(obs))
  supsim <- apply(abs(sim), 2, max)
  suppval <- sum(supsim >= supobs) / nsim
  absobs <- sum(abs(obs))
  abssim <- apply(abs(sim), 2, sum)
  abspval <- sum(abssim >= absobs) / nsim
  s2obs <- sum(obs^2)
  s2sim <- apply(sim^2, 2, sum)
  s2pval <- sum(s2sim >= s2obs) / nsim
  list(tis = fit1$tis, obs = obs, sim = sim,
       supobs = supobs, supsim = supsim, suppval = suppval,
       absobs = absobs, abssim = abssim, abspval = abspval,
       s2obs = s2obs, s2sim = s2sim, s2pval = s2pval)
}

## gofTestLinT <- function(fit1, fit2, nsim) {
## #### fit1 is the H0 model, reduced (constant coef) model
## #### fit2 is the full model
##   infl1 <- fit1$inflTheta[1,]
##   n <- length(infl1)
##   infl1 <- outer(infl1, fit1$tis)
##   p <- ncol(fit2$alpha)
##   infl2 <- sapply(fit2$inflAlpha, function(x) x[p,])
##   infl <- infl2 - infl1
##   sig <- apply(infl, 2, sd, na.rm = TRUE) / sqrt(n)
##   obs <- (fit2$alpha[,p] - fit1$theta[1] * fit1$tis) / sig
##   sim <- .Call("gofSample_rap", infl, sig, as.integer(nsim), PACKAGE="tpr")
##   supobs <- max(abs(obs))
##   supsim <- apply(abs(sim), 2, max)
##   suppval <- sum(supsim >= supobs) / nsim
##   absobs <- sum(abs(obs))
##   abssim <- apply(abs(sim), 2, sum)
##   abspval <- sum(abssim >= absobs) / nsim
##   s2obs <- sum(obs^2)
##   s2sim <- apply(sim^2, 2, sum)
##   s2pval <- sum(s2sim >= s2obs) / nsim
##   list(tis = fit1$tis, obs = obs, sim = sim,
##        supobs = supobs, supsim = supsim, suppval = suppval,
##        absobs = absobs, abssim = abssim, abspval = abspval,
##        s2obs = s2obs, s2sim = s2sim, s2pval = s2pval)
## }

plotGof <- function(gof, nuse=50, ...) {
  ran <- range(unlist(gof[2:3]))
  plot(gof$tis, gof$obs, type="s", ylim=ran, ...)
  apply(gof$sim[,1:nuse], 2,
        function(x) lines(gof$tis, x, type="s", col="gray"))
  lines(gof$tis, gof$obs, type="s")
}
