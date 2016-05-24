#############################################
#This function computes                     #
#the log  of the prior distribution of theta#
#############################################

logprior <- function(parms, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, indQ, aQ, bQ, 
                    indG, aG, bG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc){

  
  lpBeta <-   switch(indBeta, dunif(parms[1], aBeta, bBeta, log=TRUE), dgamma(parms[1], aBeta, bBeta, log=TRUE), 
                              dexp(parms[1], aBeta, log=TRUE), dnorm(parms[1], aBeta, bBeta, log=TRUE), 
                              dt(parms[1], aBeta, bBeta, log=TRUE), dweibull(parms[1], aBeta, bBeta, log=TRUE),
                              dchisq(parms[1], aBeta, bBeta, log=TRUE), dcauchy(parms[1], aBeta, bBeta, log=TRUE), 
                              dlnorm(parms[1], aBeta, bBeta,log=TRUE))


  lpQ <-   switch(indQ, dunif(parms[2], aQ, bQ, log=TRUE), dgamma(parms[2], aQ, bQ, log=TRUE), 
                              dexp(parms[2], aQ, log=TRUE), dnorm(parms[2], aQ, bQ, log=TRUE), 
                              dt(parms[2], aQ, bQ, log=TRUE), dweibull(parms[2], aQ, bQ, log=TRUE),
                              dchisq(parms[2], aQ, bQ, log=TRUE), dcauchy(parms[2], aQ, bQ, log=TRUE), 
                              dlnorm(parms[2], aQ, bQ, log=TRUE))


  lpG <-   switch(indG, dunif(parms[3], aG, bG, log=TRUE), dgamma(parms[3], aG, bG, log=TRUE), 
                              dexp(parms[3], aG, log=TRUE), dnorm(parms[3], aG, bG, log=TRUE), 
                              dt(parms[3], aG, bG, log=TRUE), dweibull(parms[3], aG, bG, log=TRUE),
                              dchisq(parms[3], aG, bG, log=TRUE), dcauchy(parms[3], aG, bG, log=TRUE), 
                              dlnorm(parms[3], aG, bG, log=TRUE))


  if(indep){
       lptauN <- tauN_sh * log(tauN_sc) - lgamma(tauN_sh) - (tauN_sh + 1) * log(parms[4]) - (tauN_sc/parms[4])
       lptauF <- tauF_sh * log(tauF_sc) - lgamma(tauF_sh) - (tauF_sh + 1) * log(parms[5]) - (tauF_sc/parms[5])
       lpvar <- lptauN + lptauF
  }
  else{
     W <- matrix(c(parms[4],parms[6],parms[6],parms[5]), 2, 2)
     lpvar <- logdiwish(W, v, S)
   }
  lp <- lpBeta + lpQ + lpG + lpvar

  return(lp)
}
