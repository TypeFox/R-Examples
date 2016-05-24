################################################################
#This function contains the Metropolis within Gibbs algorithm. # 
#The transformation if the parameters used here is the same    #
#as the one used in the BCLT algorithm.                        #
################################################################

MCMCB2zm <- function (NUpd, burnin, lag, initial, S, v, 
                        tauN_sh, tauN_sc, tauF_sh, tauF_sc, VN, 
                        VF, times, Y, indep, Sigma.Cand, 
                         m, cred,indBeta, aBeta, bBeta,
                         indQ, aQ, bQ, indG, aG, bG){

   n <- length(times)  
   if(is.null(Sigma.Cand) | is.null(initial)){
      cat("1 - Estimating initial values for Beta, Q and G... \n\n")
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
         stop("Sorry, the initial values could not be estimated. Please provide 
             initial values for Beta, Q and G.")
      }

      cat("\n\n  OK!\n\n")
      flush.console()

      cat("2 - Estimating the posterior covariance matrix for the proposal distribution...")
      flush.console()

      try(sm <- hessian(func = minus_logpost_transf, x = nlmobj$par,
         method="Richardson", method.args=list(eps=1e-5, d=0.01,
         zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2), 
         indep=indep, Y=Y, times=times, VN=VN,  VF=VF, n=n,  indBeta=indBeta, aBeta=aBeta,  
         bBeta=bBeta, rangeBeta = rangeBeta, indQ=indQ, aQ=aQ, bQ=bQ, 
         rangeQ = rangeQ, indG=indG, aG=aG, bG=bG, rangeG = rangeG, 
         S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, tauF_sh=tauF_sh, tauF_sc=tauF_sc),silent=TRUE)

      if(length(which(is.nan(sm)==TRUE))> 0){
         stop("Sorry, it was not possible to estimate the covariance matrix. 
             Please provide the covariance matrix for the proposal distribution.")
      }

      try(covMat <- solve(sm),silent=TRUE)

      if(is.null(covMat) | length(which(diag(covMat)<0))>0 | !is.pos.def(covMat)){
         stop("Sorry, it was not possible to estimate the covariance matrix. 
             Please provide the covariance matrix for the proposal distribution.")
      }

      cat("  OK!\n\n")
      flush.console()

      initial_transf <- nlmobj$par[1:3]
      Sigma.Cand_transf <- covMat[1:3,1:3]
     
      cat("3 - Running Metropolis within Gibbs...\n\n")
      flush.console()
   }
   else{
      Betas <- switch(indBeta, runif(m, aBeta, bBeta), rgamma(m, aBeta, bBeta), rexp(m, aBeta),
                rnorm(m, aBeta, bBeta), rt(m, aBeta, bBeta), rweibull(m, aBeta, bBeta),
                rchisq(m, aBeta, bBeta), rcauchy(m, aBeta, bBeta), dlnorm(m, aBeta, bBeta))

      Qs <- switch(indQ, runif(m, aQ, bQ), rgamma(m, aQ, bQ), rexp(m, aQ),
                rnorm(m, aQ, bQ), rt(m, aQ, bQ), rweibull(m, aQ, bQ),
                rchisq(m, aQ, bQ), rcauchy(m, aQ, bQ), dlnorm(m, aQ, bQ))

      Gs <- switch(indG, runif(m, aG, bG), rgamma(m, aG, bG), rexp(m, aG),
                rnorm(m, aG, bG), rt(m, aG, bG), rweibull(m, aG, bG),
                rchisq(m, aG, bG), rcauchy(m, aG, bG), dlnorm(m, aG, bG))

      rangeBeta <- range(Betas)
      rangeQ <- range(Qs)
      rangeG <- range(Gs)

      if(initial[1] < rangeBeta[1] | initial[1] > rangeBeta[2]){
         stop("Initial value for Beta is outside from its prior distribution support.")
      }  
      if(initial[2] < rangeQ[1] | initial[2] > rangeQ[2]){
         stop("Initial value for Q is outside from its prior distribution support.")
      }  
      if(initial[3] < rangeG[1] | initial[3] > rangeG[2]){
         stop("Initial value for G is outside from its prior distribution support.")
      }  

      initial_transf <- numeric()
      initial_transf[1] <- log((initial[1] - rangeBeta[1])/(rangeBeta[2] - initial[1]))
      initial_transf[2] <- log((initial[2] - rangeQ[1])/(rangeQ[2] - initial[2]))
      initial_transf[3] <- log((initial[3] - rangeG[1])/(rangeG[2] - initial[3]))
      
      saux <- rmvnorm(m, initial[1:3], Sigma.Cand)
      pos <- which(saux[,1] < rangeBeta[1] | saux[,1] > rangeBeta[2] |
                   saux[,2] < rangeQ[1] | saux[,2] > rangeQ[2] |
                   saux[,3] < rangeG[1] | saux[,3] > rangeG[2])

      saux <- saux[-pos,]
      saux_transf <- matrix(0,nrow(saux),3) 
      saux_transf[,1] <- log((saux[,1] - rangeBeta[1])/(rangeBeta[2] - saux[,1]))
      saux_transf[,2] <- log((saux[,2] - rangeQ[1])/(rangeQ[2] - saux[,2]))
      saux_transf[,3] <- log((saux[,3] - rangeG[1])/(rangeG[2] - saux[,3]))

      Sigma.Cand_transf <- cov(saux_transf)
      
     cat("1 - Running Metropolis withing Gibbs...\n\n")
     flush.console()
   }
 
   Betas <- numeric()
   Qs <- numeric()
   Gs <- numeric()
   TauNs <- numeric()
   TauFs <- numeric()
   TauNFs <- numeric()

   Betas[1] <- (rangeBeta[1] + rangeBeta[2]*exp(initial_transf[1]))/(1+ exp(initial_transf[1]))
   Qs[1] <- (rangeQ[1] + rangeQ[2]*exp(initial_transf[2]))/(1+ exp(initial_transf[2]))
   Gs[1] <- (rangeG[1] + rangeG[2]*exp(initial_transf[3]))/(1+ exp(initial_transf[3]))

   parms_transf_current <- numeric()
   parms_transf_current[1:3] <- initial_transf

   pb <- txtProgressBar(min = 1, max = NUpd, style = 3)

   for(i in 2:NUpd){
      Ytild <- log(compute_CNCF(Betas[i-1], Qs[i-1], Gs[i-1], VN, VF, times)) 
      if(indep){
         aNstar <- tauN_sh + n/2
         bNstar <- tauN_sc + 0.5*sum((Y[,1]-Ytild[,1])^2)

         aFstar <- tauF_sh + n/2
         bFstar <- tauF_sc + 0.5*sum((Y[,2]-Ytild[,2])^2)

         TauNs[i] <- 1/rgamma(1, aNstar, bNstar)
         TauFs[i] <- 1/rgamma(1, aFstar, bFstar)

         tauN_transf <- log(TauNs[i])
         tauF_transf <- log(TauFs[i])

         parms_transf_current[4:5] <- c(tauN_transf, tauF_transf)
      }
      else{
         S1 <- S + crossprod(Y-Ytild)      
         v1 <- v + n
         newSigma <- riwish(1, v1, S1)

         TauNs[i] <- newSigma[1]
         TauFs[i] <- newSigma[2]
         TauNFs[i] <- newSigma[3]

         tauN_transf <-  log(TauNs[i])
         tauF_transf <-  log(TauFs[i])
         tauNF_transf <- log((TauNFs[i] + sqrt(TauNs[i]*TauFs[i]))/(sqrt(TauNs[i]*TauFs[i]) - TauNFs[i]))

         parms_transf_current[4:6] <- c(tauN_transf, tauF_transf, tauNF_transf)
      }

      cand_transf <- rmvnorm(1, parms_transf_current[1:3], Sigma.Cand_transf)
      parms_transf_candidate <- c(cand_transf,  tauN_transf,  tauF_transf)

      if(!indep){
          parms_transf_candidate <- c(parms_transf_candidate, tauNF_transf) 
      }
     
      lp_cand <- logpost_transf(parms_transf_candidate, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, rangeBeta, indQ, aQ, bQ, rangeQ, 
                    indG, aG, bG, rangeG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc)

      lp_current <- logpost_transf(parms_transf_current, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, rangeBeta, indQ, aQ, bQ, rangeQ, 
                    indG, aG, bG, rangeG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc)


      ratio <- exp(lp_cand - lp_current)
      if(ratio > 1){
         parms_transf_current <- parms_transf_candidate
         Betas[i] <- (rangeBeta[1] + rangeBeta[2]*exp(parms_transf_current[1]))/(1+ exp(parms_transf_current[1]))
         Qs[i] <- (rangeQ[1] + rangeQ[2]*exp(parms_transf_current[2]))/(1+ exp(parms_transf_current[2]))
         Gs[i] <- (rangeG[1] + rangeG[2]*exp(parms_transf_current[3]))/(1+ exp(parms_transf_current[3]))
      }
      else{
         u <- runif(1)
         if(u < ratio){
            parms_transf_current <- parms_transf_candidate
            Betas[i] <- (rangeBeta[1] + rangeBeta[2]*exp(parms_transf_current[1]))/(1+ exp(parms_transf_current[1]))
            Qs[i] <- (rangeQ[1] + rangeQ[2]*exp(parms_transf_current[2]))/(1+ exp(parms_transf_current[2]))
            Gs[i] <- (rangeG[1] + rangeG[2]*exp(parms_transf_current[3]))/(1+ exp(parms_transf_current[3]))
         }
         else{
            Betas[i] <- Betas[i-1]
            Qs[i] <- Qs[i-1]
            Gs[i] <- Gs[i-1]
         }
      }
      setTxtProgressBar(pb, i)
   }   
   close(pb)

   cont <- length(table(Betas))
   AR = round(cont/NUpd,5)

   seq_index <- seq(burnin, NUpd, by = lag)

   ESS <- c(effectiveSize(Betas),effectiveSize(Qs), effectiveSize(Gs))

   #computing DIC
   PARAMS <- cbind(Betas[seq_index], Qs[seq_index], Gs[seq_index], TauNs[seq_index], TauFs[seq_index])
   parms <- apply(PARAMS, 2, mean)

   if(!indep) {
      PARAMS <- cbind(PARAMS, TauNFs[seq_index])
      parms <- c(parms,mean(TauNFs[seq_index]))
   }


   l <- apply(PARAMS, 1, loglik, indep = indep, Y = Y, times = times, VN = VN, VF = VF, n = n)

   Dbar <- -2*mean(l)
   Dthetabar <- -2*loglik(parms, indep, Y, times, VN, VF, n)

   pD = Dbar - Dthetabar
   DIC = pD + Dbar

   if(indep){
      r <- list(Beta=Betas, Q=Qs, G=Gs, tauN=TauNs,
          tauF=TauFs, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, 
          indep=indep, times=times, cred=cred, ESS=ESS, AR=AR,
          NUpd=NUpd, burnin=burnin, lag=lag, VN = VN, VF = VF)
    }
   else {
      r <- list(Beta=Betas, Q=Qs, G=Gs, tauN=TauNs,
          tauF=TauFs, tauNF = TauNFs, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, 
          indep=indep, times=times, cred=cred, ESS=ESS, AR=AR,
          NUpd=NUpd, burnin=burnin, lag=lag, VN = VN, VF = VF)
    }

    
   cat("\n  DONE!\n\n")
   flush.console()

   attr(r, "class") <- "mcmc"
 
   return(r)
}
