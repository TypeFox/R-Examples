#############################################################################################################
#This function is the log of a proportional function of the posterior                                       #
#distribution for a multivariate transformation of the parameters.                                          #
#This is done, so the range of the each variable is (-Inf,Inf), therefore                                   #
#the Bayesian Central Limit Theorem works properly. The most important in this function                     #
#is to know which transformation I used, and consequently how I computed the jacobian.                      #
#Suppose someone chooses priors for beta, Q and G such that the range of each parameter                     #
#is respectively rangeBeta, rangeQ and rangeG. Then 'parms_transf' is a vector such that:                   #
#                                                                                                           # 
#parms_transf[1] = log((Beta - rangeBeta[1])/(rangeBeta[2]-Beta))                                           #
#parms_transf[2] = log((Q - rangeQ[1])/(rangeQ[2]-Q))                                                       #
#parms_transf[3] = log((G - rangeG[1])/(rangeG[2]-G))                                                       #
#parms_transf[4] = log(tauN)                                                                                #
#parms_transf[5] = log(tauF)                                                                                #
#parms_transf[6] = log((tauNF + sqrt(tauN*tauF))/(sqrt(tauN*tauF)-tauNF)) (if dependent model is chosen)    #
#                                                                                                           #
#############################################################################################################


logpost_transf <- function(parms_transf, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, rangeBeta, indQ, aQ, bQ, rangeQ, 
                    indG, aG, bG, rangeG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc){

   Beta <- (rangeBeta[1] + rangeBeta[2]*exp(parms_transf[1]))/(1+ exp(parms_transf[1]))
   Q <- (rangeQ[1] + rangeQ[2]*exp(parms_transf[2]))/(1+ exp(parms_transf[2])) 
   G <- (rangeG[1] + rangeG[2]*exp(parms_transf[3]))/(1+ exp(parms_transf[3]))
   tauN <- exp(parms_transf[4])
   tauF <- exp(parms_transf[5])

   if(indep){
      parms <- c(Beta, Q, G, tauN, tauF)
   }
   else{
      tauNF <- exp(0.5*(parms_transf[4]+parms_transf[5]))*(exp(parms_transf[6])-1)/(1+exp(parms_transf[6]))
      parms <- c(Beta, Q, G, tauN, tauF, tauNF)
   }

   logjacBeta <- log((rangeBeta[2]-rangeBeta[1])*exp(parms_transf[1])/(1+exp(parms_transf[1]))^2)
   logjacQ <- log((rangeQ[2]-rangeQ[1])*exp(parms_transf[2])/(1+exp(parms_transf[2]))^2)
   logjacG <- log((rangeG[2]-rangeG[1])*exp(parms_transf[3])/(1+exp(parms_transf[3]))^2)
   logjacobian <- logjacBeta + logjacQ + logjacG + parms_transf[4] + parms_transf[5]

   if(!indep){
      logjacobian <- logjacobian + 0.5*(parms_transf[4] + parms_transf[5]) + parms_transf[6] + log(2) - 2*log(1+exp(parms_transf[6]))
   }

   lp <- logpost(parms, indep, Y, times, VN, VF, n, 
             indBeta, aBeta, bBeta, indQ, aQ, bQ, 
             indG, aG, bG, S, v, tauN_sh, tauN_sc, tauF_sh, tauF_sc) + logjacobian

   return(lp)
   }
