###############################################
#This function computes the minus of the log  #
#posterior of the transformed parameters      #
###############################################

minus_logpost_transf <- function(parms_transf, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, rangeBeta, indQ, aQ, bQ, rangeQ, 
                    indG, aG, bG, rangeG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc){

   mlp <- -logpost_transf(parms_transf, indep, Y, times, VN, VF, n, 
                    indBeta, aBeta, bBeta, rangeBeta, indQ, aQ, bQ, rangeQ, 
                    indG, aG, bG, rangeG, S, v, 
                    tauN_sh, tauN_sc, tauF_sh, tauF_sc)

   return(mlp)
   }
