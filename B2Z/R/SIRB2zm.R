#########################################################
#This function implements the SIR Algorithm             #
#indBeta, indQ and indG are numbers that indicates which#
#prior distribution the user has chosen.                #
#########################################################

SIRB2zm <- function (m, Y, times, VN, VF, indBeta, aBeta, bBeta, 
                     indQ, aQ, bQ, indG, aG, bG, v, S, tauN_sh, 
                     tauN_sc, tauF_sh, tauF_sc, indep, cred){

    n <- length(times)

    cat("1 - Sampling from the prior distribution of theta...")
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
 
    cat("  OK!\n\n")
    flush.console()
 
   
    l <- apply_pb(PARAMS, 1, loglik, indep = indep, Y = Y, times = times, VN = VN, VF = VF, n = n)

    C <- 700 - (max(range(l)))
    weights <- l + C
    weights <- exp(weights)
    weights <- weights/sum(weights)

    draw_index <- sample(1:m, m, replace=TRUE, prob=weights)
   
    Betaout <- Betas[draw_index]
    Qout <- Qs[draw_index]
    Gout <- Gs[draw_index]
    TauNout  <- TauNs[draw_index]
    TauFout  <- TauFs[draw_index]
    
    if(!indep){TauNFout <- TauNFs[draw_index]}
    else{TauNFout <- 0}

    #Proportion of diferent values
    prop <- round(length(table(Betaout))/m,6)

    #Computing DIC 
    Dbar <- -2*mean(l[draw_index])
    parms <- c(mean(Betaout),mean(Qout),mean(Gout), mean(TauNout), mean(TauFout), mean(TauNFout))

    Dthetabar <- -2*loglik(parms, indep, Y, times, VN, VF, n)

    pD = Dbar - Dthetabar
    DIC = pD + Dbar

   #Computing ESS
    ESS = m/(1 + var(weights))

   if(indep)
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauF=TauFout, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, 
          ESS=ESS, prop=prop, indep=indep, times=times, 
          weights=weights, maxw= max(weights), cred=cred, 
          VN = VN, VF = VF)
     }
   else 
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = TauNFout, tauF=TauFout, Y=Y, DIC=DIC, 
          pD=pD, Dbar=Dbar, ESS=ESS, prop=prop, indep=indep, 
          times=times, weights=weights, maxw= max(weights),
          cred=cred, VN = VN, VF = VF)
     }

   cat("\nDONE!\n\n")
   flush.console()

   attr(r, "class") <- "sir"
   return(r)

}
