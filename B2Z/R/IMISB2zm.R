##############################################
#This function implements the IMIS algorithm #
##############################################

IMISB2zm <- function (N0, B, M, it.max, S, v, 
              tauN_sh, tauN_sc, tauF_sh, tauF_sc,
              VN, VF, times, Y, indep, cred,
              indBeta, aBeta, bBeta,
              indQ, aQ, bQ, indG, aG, bG){

   n <- length(times)

   cat("1 - Initial Stage...\n\n")
   flush.console()

   Betas <- switch(indBeta, runif(N0+2, aBeta, bBeta), rgamma(N0+2, aBeta, bBeta), rexp(N0+2, aBeta),
                rnorm(N0+2, aBeta, bBeta), rt(N0+2, aBeta, bBeta), rweibull(N0+2, aBeta, bBeta),
                rchisq(N0+2, aBeta, bBeta), rcauchy(N0+2, aBeta, bBeta), dlnorm(N0+2, aBeta, bBeta))

   Qs <- switch(indQ, runif(N0+2, aQ, bQ), rgamma(N0+2, aQ, bQ), rexp(N0+2, aQ),
                rnorm(N0+2, aQ, bQ), rt(N0+2, aQ, bQ), rweibull(N0+2, aQ, bQ),
                rchisq(N0+2, aQ, bQ), rcauchy(N0+2, aQ, bQ), dlnorm(N0+2, aQ, bQ))

   Gs <- switch(indG, runif(N0+2, aG, bG), rgamma(N0+2, aG, bG), rexp(N0+2, aG),
                rnorm(N0+2, aG, bG), rt(N0+2, aG, bG), rweibull(N0+2, aG, bG),
                rchisq(N0+2, aG, bG), rcauchy(N0+2, aG, bG), dlnorm(N0+2, aG, bG))


   if (length(which(Betas < 0)) > 0) {
      stop("Negative values for Beta were generated from its prior distribution. Try another prior distribution.")
   }
   if (length(which(Qs < 0)) > 0) {
      stop("Negative values for Q were generated from its prior distribution. Try another prior distribution.")
   }
   if (length(which(Gs < 0)) > 0) {
      stop("Negative values for G were generated from its prior distribution. Try another prior distribution.")
   }

   rangeBeta <- range(Betas)
   rangeQ <- range(Qs)
   rangeG <- range(Gs)

   posmax <- which.max(Betas)
   posmin <- which.min(Betas)
   Betas <- Betas[-c(posmax,posmin)]
   
   posmax <- which.max(Qs)
   posmin <- which.min(Qs)
   Qs <- Qs[-c(posmax,posmin)]
  
   posmax <- which.max(Gs)
   posmin <- which.min(Gs)
   Gs <- Gs[-c(posmax,posmin)]


   if(indep){
      TauNs <- 1/rgamma(N0, tauN_sh, tauN_sc)
      TauFs <- 1/rgamma(N0, tauF_sh, tauF_sc)
      PARAMS <- cbind(Betas, Qs, Gs, TauNs, TauFs)
   }
   else{
      Taus <- sapply(1:N0, riwish, v = v, S = S)
      TauNs <- Taus[1,]
      TauFs <- Taus[2,]
      TauNFs <- Taus[3,]
      PARAMS <- cbind(Betas, Qs, Gs, TauNs, TauFs, TauNFs)
   }
 
   L <- apply_pb(PARAMS, 1, loglik, indep = indep, Y = Y, times = times, VN = VN, VF = VF, n = n)

   C <- 700 - (max(range(L)))
   w <- L + C
   w <- exp(w)
   w <- w/sum(w)

   params_transf <- matrix(0,nrow(PARAMS),ncol(PARAMS))
   params_transf[,1] <- log((PARAMS[,1] - rangeBeta[1])/(rangeBeta[2] - PARAMS[,1]))
   params_transf[,2] <- log((PARAMS[,2] - rangeQ[1])/(rangeQ[2] - PARAMS[,2]))
   params_transf[,3] <- log((PARAMS[,3] - rangeG[1])/(rangeG[2] - PARAMS[,3]))
   params_transf[,4:5] <- log(PARAMS[,4:5])

   if(!indep){
      params_transf[,6] <- log((PARAMS[,6] + sqrt(PARAMS[,4]*PARAMS[,5]))/(sqrt(PARAMS[,4]*PARAMS[,5])-PARAMS[,6]))
   } 

   ptheta <- exp(apply_pb( params_transf, 1, logprior_transf, indep = indep, Y = Y, times = times, VN = VN, VF =VF, 
                  n = n, indBeta = indBeta, aBeta = aBeta, bBeta = bBeta, rangeBeta = rangeBeta, 
                  indQ = indQ, aQ = aQ, bQ = bQ, rangeQ = rangeQ, indG = indG, aG = aG, bG = bG, 
                  rangeG = rangeG, S = S, v = v, tauN_sh = tauN_sh, tauN_sc = tauN_sc, 
                  tauF_sh = tauF_sh, tauF_sc = tauF_sc))

   cat("\nOK!\n\n")
   flush.console()

 
  #Algorithm Progression 
   plot(0,0, xlab = "Iteration", main = "IMIS Algorithm Progression", ylab ="Expected fraction of unique points", 
         pch = 19, ylim=c(0,1), xlim = c(0, 1.2*it.max), axes=F)

   box()
   abline(h=1-1/exp(1), lwd=2, col = "red")
   axis(2,at=(1-1/exp(1)), labels = "0.632")
   axis(1, at=seq(it.max,1,l=it.max/3), labels = c("Max",
   round(seq(it.max,1,l=it.max/3)[-1],1)))
   abline(v=it.max, lwd=2, col = "blue")

   #Importance Sampling Stage
   cat("\n2 - Importance Sampling Stage...  ")
   flush.console()
   dimens <- ifelse(indep, 5, 6)

   expfrac <- numeric()
   expfrac[1] <- 0
   k <- 0
   theta <- list()  
   SigmaK <- list()

   thetas <- params_transf
   
   Qi <- 0

   while(Qi < (1-exp(-1)) & k <= (it.max-1) ){
      #a
      k <- k + 1
      Nk <- N0 + B*k
      ind.max   <- which.max(w)
      theta[[k]] <- thetas[ind.max,]
      inv.cov <- solve(cov(thetas))
      distances <- apply(thetas, 1, mahalanobis_mod, d = dimens, center = theta[[k]], invcov = inv.cov)
      o <- order(distances)[1:B]
      Mat <- thetas[o,]
      wt <- (w[o] + 1/Nk)/2
      wt <- wt/sum(wt)
      SigmaK[[k]] <- cov.wt(Mat, wt = wt)$cov

      #b
      newinput <- rmvnorm(B, theta[[k]], SigmaK[[k]])
      thetas <- rbind(thetas, newinput)
  

      #c

      newBetas <-  (rangeBeta[1] + rangeBeta[2]*exp(newinput[,1]) )/ (1 + exp(newinput[,1]))
      newQs <- (rangeQ[1] + rangeQ[2]*exp(newinput[,2]) )/ (1 + exp(newinput[,2]))
      newGs <- (rangeG[1] + rangeG[2]*exp(newinput[,3]) )/ (1 + exp(newinput[,3]))
      newtauNs <- exp(newinput[,4])
      newtauFs <- exp(newinput[,5])

      PARAMS <- cbind(newBetas, newQs, newGs, newtauNs, newtauFs)

      if(!indep){
         newtauNFs <- exp(0.5*(newinput[,4]+newinput[,5]))*(exp(newinput[,6])-1)/(1+exp(newinput[,6]))
         PARAMS <- cbind(PARAMS, newtauNFs)
      }


      newL <- apply(PARAMS, 1, loglik, indep = indep, Y = Y, times = times, VN = VN, VF = VF, n = n)
      L <- c(L, newL)


      newptheta <- exp(apply(newinput, 1, logprior_transf, indep = indep, Y = Y, times = times, VN = VN, VF =VF, 
                  n = n, indBeta = indBeta, aBeta = aBeta, bBeta = bBeta, rangeBeta = rangeBeta, 
                  indQ = indQ, aQ = aQ, bQ = bQ, rangeQ = rangeQ, indG = indG, aG = aG, bG = bG, 
                  rangeG = rangeG, S = S, v = v, tauN_sh = tauN_sh, tauN_sc = tauN_sc, 
                  tauF_sh = tauF_sh, tauF_sc = tauF_sc))


      ptheta <- c(ptheta, newptheta)
        if(k==1) 
          {
          invsigma <- solve(SigmaK[[1]])
          sumHs <- apply(thetas, 1, dmnorm_mod, mean=theta[[1]], d=dimens, varcov=SigmaK[[1]], invvarcov=invsigma)
          }
        else
          {
          sumHsnew <- rep(0,B)
          for(s in 1:(k-1))
             {
             invsigma <- solve(SigmaK[[s]])
             sumHsnew <- sumHsnew + apply(newinput,1, dmnorm_mod, mean=theta[[s]],d=dimens,varcov=SigmaK[[s]], invvarcov=invsigma)
             }
           
  
          sumHs <- c(sumHs, sumHsnew)
          invsigma <- solve(SigmaK[[k]])
          sumHs <- sumHs + apply(thetas, 1, dmnorm_mod, mean=theta[[k]],d=dimens, varcov=SigmaK[[k]], invvarcov=invsigma)
          }  

        qtheta <- (N0/Nk)*ptheta + (B/Nk)*sumHs
        C <- 700 - (max(range(L)))
        w <- exp(L + C)*ptheta/qtheta
        w <- w/sum(w)
        Qi <- sum(1-(1-w)^M)/M
        expfrac[k+1] <- Qi
   
        segments((k-1), expfrac[k], k, expfrac[k+1])
        points(k, expfrac[k+1],pch=19)
        }

   cat("   OK!\n\n")
   flush.console()
   
   cat("3 - Resample Stage...")
   flush.console()

   V.hat <- sum((Nk*w-1)^2)/Nk
   U.hat <-    -log(prod(w^(w/log(Nk))))
   Q.hat <- sum(1-(1-w)^M)
   ESS <- 1/sum(w^2)
   maxw <- max(w)
   options(warn=-1)
   draw_index <- sample(1:Nk, M, replace=TRUE, prob=w)
  
   options(warn=0)
   if(k==(it.max+1)){
   warning("Expected fraction of unique points < (1-1/e)")}

   Betaout <-  (rangeBeta[1] + rangeBeta[2]*exp(thetas[draw_index,1]) )/ (1 + exp(thetas[draw_index,1]))
   Qout <- (rangeQ[1] + rangeQ[2]*exp(thetas[draw_index,2]) )/ (1 + exp(thetas[draw_index,2]))
   Gout <- (rangeG[1] + rangeG[2]*exp(thetas[draw_index,3]) )/ (1 + exp(thetas[draw_index,3]))
   TauNout <- exp(thetas[draw_index,4])
   TauFout <- exp(thetas[draw_index,5])

   parms <- c(mean(Betaout),mean(Qout),mean(Gout), mean(TauNout), mean(TauFout))

   if(!indep){
      TauNFout <- exp(0.5*(thetas[draw_index,4]+thetas[draw_index,5]))*(exp(thetas[draw_index,6])-1)/(1+exp(thetas[draw_index,6]))
      parms <- c(parms, mean(TauNFout))
      PARAMS <- cbind(PARAMS, newtauNFs)
   }

  
   wout <- w[draw_index]

   Dbar <- -2*mean((L[draw_index]))
   Dthetabar <- -2*loglik(parms, indep, Y, times, VN, VF, n)

   pD = Dbar - Dthetabar
   DIC = pD + Dbar

   if(indep){
      r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = NULL, tauF=TauFout, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, ESS=ESS, Qi=Qi, 
          indep=indep, times=times, cred=cred, 
          expfrac=expfrac, V.hat=V.hat, U.hat=U.hat, Q.hat=Q.hat, 
          maxw=maxw, w=wout, VN = VN, VF = VF)
   }
   else {
      r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = TauNFout, tauF=TauFout, Y=Y,
          DIC=DIC, ESS=ESS, pD=pD, Dbar=Dbar, Qi=Qi, indep=indep,
          times=times, cred=cred, 
          expfrac=expfrac, V.hat=V.hat, U.hat=U.hat, Q.hat=Q.hat, 
          maxw=maxw, w=wout, VN = VN, VF = VF)
   }

   cat("   OK!\n")
   flush.console()
   cat("\nDONE!\n")
   flush.console()

   attr(r, "class") <- "imis"
 
   return(r)
   }

