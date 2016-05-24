###
### pplProbit.R
###

##-----------------------------------------------------------------------------
## Find probit model minimizing Posterior Predictive Loss
pplPM <- function(tauseq=exp(seq(log(.01), log(2), length=20)),
                  kPen=1,
                  y,
                  x,
                  xadj,
                  niter=10^4,
                  thinning=1,
                  burnin=round(niter/10),
                  priorCoef,
                  priorDelta,
                  priorVar,
                  initSearch='greedy',
                  mc.cores=1) {
  if (missing(priorDelta)) {
    priorDelta <- new("msPriorSpec",
                      priorType='modelIndicator',
                      priorDistr='uniform',
                      priorPars=double(0))
  }
  if (missing(priorVar)) {
    priorVar <- new("msPriorSpec",
                    priorType='nuisancePars',
                    priorDistr='invgamma',
                    priorPars=c(alpha=.01, lambda=.01))
  }
  if (priorCoef@priorDistr == 'pMOM') {
    nlpfit <- function(tau) {
      pr <- priorCoef
      pr@priorPars['tau'] <- tau
      fit1 <- pmomPM(y=y,
                     x=x,
                     xadj=xadj,
                     niter=niter,
                     thinning=thinning,
                     burnin=burnin,
                     priorCoef=pr,
                     priorDelta=priorDelta,
                     initSearch=initSearch,
                     verbose=FALSE)
      p1 <- pplProbit(fit1, x=x, xadj=xadj, y=y, kPen=kPen)
      list(fit1=fit1,
           ppl=c(p1$d, p1$p, p1$g, p1$msize))
    }
  } else if (priorCoef@priorDistr == 'peMOM') {
    nlpfit <- function(tau) {
      pr <- priorCoef
      pr@priorPars['tau'] <- tau
      fit1 <- emomPM(y=y,
                     x=x,
                     xadj=xadj,
                     niter=niter,
                     thinning=thinning,
                     burnin=burnin,
                     priorCoef=pr,
                     priorDelta=priorDelta,
                     initSearch=initSearch,
                     verbose=FALSE)
      p1 <- pplProbit(fit1, x=x, xadj=xadj, y=y, kPen=kPen)
      list(fit1=fit1,
           ppl=c(p1$d, p1$p, p1$g, p1$msize))
    }
  } else {
    stop('prior on coefficients not recognized. Currently only pMOM and peMOM are implemented')
  }

  nlpseq <- if ("parallel" %in% loadedNamespaces())  {
                       parallel::mclapply(as.list(tauseq),
                       nlpfit,
                       mc.cores=mc.cores,
                       mc.preschedule=FALSE)
            }
            else {
              lapply(as.list(tauseq), nlpfit)
            }

  PPL <- data.frame(tauseq, do.call(rbind, lapply(nlpseq, '[[', 'ppl')))
  colnames(PPL) <- c('tau', 'PPL', 'G', 'P', 'msize')
  #require(mgcv)
  gam1 <- gam(PPL ~ s(tau, k=min(10, length(tauseq))),
              family=Gamma(link="log"),
              data=PPL,
              sp= -1)
  PPL$sPPL <- exp(predict(gam1))
  ans <- list(nlpseq[[which.min(PPL$sPPL)]]$fit1,
              PPL,
              PPL$tau[which.min(PPL$PPL)])
  names(ans) <- c('optfit', 'PPL', 'tauopt')
  ans
}


##-----------------------------------------------------------------------------
## Evaluate Posterior Predictive Loss under a probit model.
## Input:
## - fit: probit model fit, e.g. as returned by pmomPM
## - x: covariates used to fit the model
## - xadj: adjustment covariates
## - y: response variables (e.g. 0/1 vector or TRUE/FALSE)
## - kPen: Loss is Dev(yp,a) + kPen*Dev(yobs,a), where yp: draw from post predictive, yobs: observed data and a is E(yp|yobs).
##         i.e. kPen is a penalty term specifying the relative importance of deviations from the observed data.
## Ouput
## - D: P + G
## - P: P_k(m) - (Penalty), i.e. sum(hm - h(mu))
## - G: G_k(m) - (Fit)
pplProbit <- function(fit, x, xadj, y, kPen=1) {
  m <- nrow(fit$postModel)
  n <- nrow(x)
  ## Generate predictive distribution -------------------------------------
  yp <- matrix(0, nrow=m, ncol=n)
  for (i in 1:m) {
    th1 <- fit$postCoef1[i, ]
    th2 <- fit$postCoef2[i, ]
    lpred <- as.vector(x %*% matrix(th1, ncol=1) + xadj %*% matrix(th2, ncol=1))
    p <- pnorm(lpred)
    yp[i, ] <- rbinom(n, 1, p)
  }
  ## Compute ppl ----------------------------------------------------------
  h <- function(z) {
    (z + 0.5)*log(z + 0.5) + (1.5-z)*log(1.5 - z)
  }
  msize <- mean(apply(fit$postModel, 1, sum)) + ncol(xadj)
  if (kPen == 'msize') {
    kPen <- msize
  }
  ## P_k(m) - (Penalty) ---------------------------------------------------
  mu    <- apply(yp, 2, mean, na.rm=TRUE)
  hi    <- h(yp)
  hm    <- apply(hi, 2, mean, na.rm=TRUE)
  P     <- sum(hm - h(mu))
  ## G_k(m) - (Fit) -------------------------------------------------------
  Gm    <- (h(mu) + kPen*h(y))/(kPen+1) - h((mu + kPen*y)/(kPen+1))
  G     <- (kPen+1)*sum(Gm, na.rm=TRUE)
  ## D_k(m) ---------------------------------------------------------------
  D <- P + G
  ## Return Output --------------------------------------------------------
  list(d=D,
       g=G,
       p=P,
       msize=msize)
}

