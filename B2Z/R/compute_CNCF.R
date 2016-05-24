################################################################
#This function computes the concentrations at the near and far #
#field at time t, given the parameters of the two-zone model:  #
#Beta, Q, G, VN, VF.                                           #
################################################################ 
 
compute_CNCF <- function(Beta, Q, G, VN, VF, t){
   lambda1 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) + sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))
   lambda2 <- 0.5*(-(Beta*VF+(Beta+Q)*VN)/(VN*VF) - sqrt(((Beta*VF+(Beta+Q)*VN)/(VN*VF))^2 - 4*(Beta*Q/(VN*VF))))
   elambda1 <- exp(lambda1*t)
   elambda2 <- exp(lambda2*t)
   CN <- G/Q + G/Beta + G*((Beta*Q + lambda2*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda1 - 
                        G*((Beta*Q + lambda1*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda2

   CF <- G/Q + G*((lambda1*VN+Beta)/Beta)*((Beta*Q + lambda2*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda1 - 
               G*((lambda2*VN+Beta)/Beta)*((Beta*Q + lambda1*VN*(Beta+Q))/(Beta*Q*VN*(lambda1-lambda2)))*elambda2 

   return(cbind(CN,CF))
}
