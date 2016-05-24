#=======================================================================================================================================#
# Bayesian Sparse Group Selection                                                                                                       #
# Date: 02/28/2013                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Determine if the selection group is added or deleted                                                                     #
#=======================================================================================================================================#


GroupMH = function(Y, X, beta.value, r, tau2, rho, sigma2, theta, eta)
{
   log.ratio.prob.comp = 0
   X = cbind(X)
   beta.value = cbind(beta.value)
   for(var.index in 1:ncol(X)){
     Ri = Y- cbind(X[, -var.index]) %*% cbind(beta.value[-var.index, ])
     ri = t(Ri)%*% cbind(X[, var.index]) * tau2[var.index] /( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
     sigma2i.star = sigma2 * tau2[var.index]/( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
     log.zi = 0.5* ( log(sigma2i.star)-log(tau2[var.index]) ) + ri^2/(2*sigma2i.star)
     if(r[var.index] == 1){
          log.prob.num = -0.5*(log(2*pi) + log(tau2[var.index])) -  ( (beta.value[var.index])^2/(2*tau2[var.index])) +
                         log(1-rho[var.index])

          log.prob.den = -0.5*(log(2*pi) + log(sigma2i.star)) - ( (beta.value[var.index]- ri)^2/(2*sigma2i.star)) +
                         log(1-rho[var.index])+ log.zi - log(rho[var.index] + (1-rho[var.index])*exp(log.zi))

     }
      else{
          log.prob.num = log(rho[var.index])
          log.prob.den = log(rho[var.index]) - log( rho[var.index] + (1-rho[var.index])*exp(log.zi))
      }

      log.ratio.prob.comp = log.ratio.prob.comp + log.prob.num - log.prob.den/ncol(X)
   }
   log.ratio.prob.comp = log(1- theta) - log(theta) - 0.5*(t(beta.value)%*%t(X)%*%X%*%beta.value-2*t(Y)%*%X%*%beta.value)/sigma2 + log.ratio.prob.comp

   if(eta==0)# 0 ==> 1
     eta.sample = ifelse(( log(runif(1)) < log.ratio.prob.comp ), 1, 0)
   else # 1 ==> 0
     eta.sample = ifelse(( log(runif(1)) < -log.ratio.prob.comp), 0, 1)

   as.integer(eta.sample)
}
GroupMH = cmpfun(GroupMH)
