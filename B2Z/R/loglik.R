##########################################################
#This function computes the log-likehood.                #
#The input parms is a vector that has 6 elements,        #
#such that: parms[1] = Beta, parms[2] = Q, parms[3] = G  #
#parms[4] = tauN, parms[5] = tauF, parms[6] = tauNF.     #
#If the independent model is considered then parms[6] is #
#always 0.                                               #
##########################################################

loglik <- function(parms, indep, Y, times, VN, VF, n){
   if(indep){
      Sigmai <- matrix(c(1/parms[4], 0, 0, 1/parms[5]), 2, 2)
   }
   else{
      Sigmai <- (1/(parms[4]*parms[5] - parms[6]^2))*matrix(c(parms[5], -parms[6], -parms[6], parms[4]), 2, 2)
   }
   Ytild <- log(compute_CNCF(parms[1], parms[2], parms[3], VN, VF, times))
   Y_Minus_Ytild <- Y - Ytild
   l <- -n*log(2*pi) + (n/2)*determinant(Sigmai)$modulus - 0.5*sum(mat_mul(Y_Minus_Ytild[,1],Y_Minus_Ytild[,2],Sigmai))
   return(l)
}