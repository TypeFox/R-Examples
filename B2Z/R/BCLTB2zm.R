##############################################################
#This function implements the Bayesian Central Limit Theorem #
#for the Bayesian -Two Zone Model (B2Z)                      #                                                            #
#indBeta, indQ and indG are numbers that indicates which     #
#prior distribution the user has chosen.                     #
##############################################################

BCLTB2zm <- function(S, v, tauN_sh, tauN_sc, tauF_sh, tauF_sc, VN, 
                     VF, times, Y, indep, m, cred, indBeta, aBeta, bBeta,
                     indQ, aQ, bQ, indG, aG, bG, size_sample){

    n <- length(times)

    cat("1 - Computing starting values... \n\n")
    flush.console()

    Betas <- switch(indBeta, runif(m, aBeta, bBeta), rgamma(m, aBeta, bBeta), rexp(m, aBeta),
                rnorm(m, aBeta, bBeta), rt(m, aBeta, bBeta), rweibull(m, aBeta, bBeta),
                rchisq(m, aBeta, bBeta), rcauchy(m, aBeta, bBeta), dlnorm(m, aBeta, bBeta))

    Qs <- switch(indQ, runif(m, aQ, bQ), rgamma(m, aQ, bQ), rexp(m, aQ),
                rnorm(m, aQ, bQ), rt(m, aQ, bQ), rweibull(m, aQ, bQ),
                rchisq(m, aQ, bQ), rcauchy(m, aQ, bQ), dlnorm(m, aQ, bQ))

    Gs <- switch(indG, runif(m, aG, bG), rgamma(m, aG, bG), rexp(m, aG),
                rnorm(m, aG, bG), rt(m, aG, bG), rweibull(m, aG, bG),
                rchisq(m, aG, bG), rcauchy(m, aG, bG), dlnorm(m, aG, bG))

    if(indep){
        TauNs <- 1/rgamma(m, tauN_sh, tauN_sc)
        TauFs <- 1/rgamma(m, tauF_sh, tauF_sc)
        PARAMS <- cbind(Betas, Qs, Gs, TauNs, TauFs)
     }
    else{
        Taus <- sapply(1:m, riwish, v = v, S = S)
        TauNs <- Taus[1,]
        TauFs <- Taus[2,]
        TauNFs <- Taus[3,]
        PARAMS <- cbind(Betas, Qs, Gs, TauNs, TauFs, TauNFs)
     }




  lp <- apply_pb(PARAMS, 1, logpost, indep=indep, Y=Y, times=times, VN=VN, VF=VF, n=n,
                 indBeta = indBeta, aBeta = aBeta, bBeta = bBeta, indQ = indQ, 
                 aQ = aQ, bQ = bQ, indG = indG, aG = aG, bG = bG, S = S, v = v, 
                 tauN_sh = tauN_sh, tauN_sc = tauN_sc, tauF_sh = tauF_sh, 
                 tauF_sc = tauF_sc)


  cat("\n2 - Finding posterior mode...") 
  flush.console()
  
  pos <- which.max(lp)
  rangeBeta <- range(PARAMS[,1])
  rangeQ <- range(PARAMS[,2])
  rangeG <- range(PARAMS[,3])

  initial <- PARAMS[pos,]
  initial_transf <- numeric()
  initial_transf[1] <- log((initial[1] - rangeBeta[1])/(rangeBeta[2] - initial[1]))
  initial_transf[2] <- log((initial[2] - rangeQ[1])/(rangeQ[2] - initial[2]))
  initial_transf[3] <- log((initial[3] - rangeG[1])/(rangeG[2] - initial[3]))
  initial_transf[4:5] <- log(initial[4:5])

  if(!indep){
     initial_transf[6] <- log((initial[6] + sqrt(initial[4]*initial[5]))/(sqrt(initial[4]*initial[5])-initial[6]))
  }

  try(nlmobj <- nlminb(start = initial_transf, objective = minus_logpost_transf, indep=indep, Y=Y,
                   times=times, VN=VN,  VF=VF, n=n,  indBeta=indBeta, aBeta=aBeta,  
                   bBeta=bBeta, rangeBeta = rangeBeta, indQ=indQ, aQ=aQ, bQ=bQ, 
                   rangeQ = rangeQ, indG=indG, aG=aG, bG=bG, rangeG = rangeG, 
                   S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, tauF_sh=tauF_sh, tauF_sc=tauF_sc), silent=TRUE)

   if(!exists("nlmobj")){
       stop("Sorry, the posterior mode could not be found. Try the other methods: IMIS, SIR or MCMC.")
   }
   cat("  OK!\n\n")
   flush.console()
   cat("3 - Estimating the posterior covariance matrix...")
   flush.console()

   try(sm <- hessian(func = minus_logpost_transf, x = nlmobj$par,
       method="Richardson", method.args=list(eps=1e-5, d=0.01,
       zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2), 
       indep=indep, Y=Y, times=times, VN=VN,  VF=VF, n=n,  indBeta=indBeta, aBeta=aBeta,  
       bBeta=bBeta, rangeBeta = rangeBeta, indQ=indQ, aQ=aQ, bQ=bQ, 
       rangeQ = rangeQ, indG=indG, aG=aG, bG=bG, rangeG = rangeG, 
       S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, tauF_sh=tauF_sh, tauF_sc=tauF_sc),silent=TRUE)

    if(length(which(is.nan(sm)==TRUE))> 0){
       stop("Sorry, it was not possible to estimate the covariance matrix. Try the other methods: IMIS, SIR or MCMC.")
    }

    try(covMat <- solve(sm),silent=TRUE)

    if(is.null(covMat) || length(which(diag(covMat)<0))>0 || !is.pos.def(covMat)){
       stop("Sorry, it was not possible to estimate the covariance matrix. Try the other methods: IMIS, SIR or MCMC.")
    }

   cat("  OK!\n\n")
   flush.console()

   cat("4 - Sampling from the posterior distribution...")

   dt <- rmvnorm(size_sample, nlmobj$par, covMat)

   Betas <-  (rangeBeta[1] + rangeBeta[2]*exp(dt[,1]) )/ (1 + exp(dt[,1]))
   Qs <- (rangeQ[1] + rangeQ[2]*exp(dt[,2]) )/ (1 + exp(dt[,2]))
   Gs <- (rangeG[1] + rangeG[2]*exp(dt[,3]) )/ (1 + exp(dt[,3]))
   TauNs <- exp(dt[,4])
   TauFs <- exp(dt[,5])

   PARAMS <- cbind(Betas, Qs, Gs, TauNs, TauFs)
   parms <- c(mean(Betas),mean(Qs),mean(Gs), mean(TauNs), mean(TauFs))
   if(!indep){
      TauNFs <- exp(0.5*(dt[,4]+dt[,5]))*(exp(dt[,6])-1)/(1+exp(dt[,6]))
      PARAMS <- cbind(PARAMS, TauNFs)
      parms <- c(parms, mean(TauNFs))
   }


    l <- apply(PARAMS, 1, loglik, indep = indep, Y = Y, times = times, VN = VN, VF = VF, n = n)

    Dbar <- -2*mean(l)
    Dthetabar <- -2*loglik(parms, indep, Y, times, VN, VF, n)

    pD = Dbar - Dthetabar
    DIC = pD + Dbar

    cat("  DONE!\n\n")
    flush.console()

    if(indep){
      r <- list(Beta=Betas, Q=Qs, G=Gs, tauN=TauNs,
          tauF=TauFs, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, 
          indep=indep, times=times, cred=cred, VN = VN, VF = VF)
    }
    else {
      r <- list(Beta=Betas, Q=Qs, G=Gs, tauN=TauNs,
          tauF=TauFs, tauNF = TauNFs, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, 
          indep=indep, times=times, cred=cred, VN = VN, VF = VF)
    }

    attr(r, "class") <- "bclt"

    return(r)
    }