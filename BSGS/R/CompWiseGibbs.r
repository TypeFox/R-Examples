#=======================================================================================================================================#
# Bayesian Variable Selection Approach with Gibbs Sampler                                                                               #
# Date: 02/28/2013                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Performa Bayesian variable selection with component-wise prior to identify the important variable                        #
#=======================================================================================================================================#

CompWiseGibbs = function(Y, X, beta.value, r, tau2, rho, sigma2, Iter)
{
  for(iter in 1:Iter){
    for(var.index in 1:ncol(X)){
      Res = Y- cbind(X[, -var.index]) %*% cbind(beta.value[-var.index, ])
      ri = t(Res)%*%cbind( X[, var.index] )* tau2[var.index] /( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
      sigma2i.star = sigma2 * tau2[var.index]/( sum(X[, var.index]^2)*tau2[var.index] + sigma2 )
      log.zi = 0.5* ( log(sigma2i.star)-log(tau2[var.index]) ) + ri^2/(2*sigma2i.star)
      prob.0 = rho[var.index]/( rho[var.index]+ (1-rho[var.index])* exp(log.zi) ) 
      r[var.index] = ifelse( runif(1)< prob.0, 0, 1)
      beta.value[var.index] = ifelse(r[var.index], rnorm(1, ri, sqrt(sigma2i.star)), 0)
    }
  }
    list(beta.sample = beta.value, r.sample = r)
}
CompWiseGibbs = cmpfun(CompWiseGibbs)
 
 
 

CompWiseGibbsSimple = function(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda, num.of.inner.iter, num.of.iteration, MCSE.Sigma2.Given)
{
	num.of.obs = length(Y)
	num.of.covariates = ncol(X)
	r.samples = sigma2.samples = beta.samples = numeric(0)

  	r.sample = r.samples = r
  	beta.sample = beta.samples = c(beta.value)
  	sigma2.sample = sigma2.samples = sigma2

  iter = 1
  start.run = proc.time()
  while( iter< num.of.iteration){
		for(num.of.inner.iter.index in 1:num.of.inner.iter){
			for(var.index in 1:ncol(X)){
				Res = Y- X[, -var.index, drop=FALSE] %*% cbind(beta.sample[-var.index])
				ri = ( t(Res)%*%X[, var.index, drop=FALSE] )* tau2 /( sum(X[, var.index]^2)*tau2 + sigma2.sample)
				sigma2i.star = sigma2.sample * tau2/( sum(X[, var.index]^2)*tau2+ sigma2.sample )
				log.zi = 0.5* ( log(sigma2i.star)-log(tau2) ) + ri^2/(2*sigma2i.star)
				prob.0 = rho/( rho+ (1-rho)* exp(log.zi) ) #(1-rho[var.index])*exp(log.zi) / (rho[var.index] + (1-rho[var.index])*exp(log.zi) )#1/(1 + 1 + (rho[var.index]/(1-rho[var.index]))* exp(log.zi)) #
				r.sample[var.index] = ifelse( (runif(1)< prob.0), 0, 1)
				beta.sample[var.index] = ifelse(r.sample[var.index], rnorm(1, ri, sqrt(sigma2i.star)), 0)
			}
		}
    Residual2 = sum((Y-X%*%matrix(beta.sample, ncol=1))^2)
    sigma2.sample = 1/rgamma(1, (length(Y)+nu)/2, (Residual2+nu*lambda)/2) 

		
	beta.samples = cbind(beta.samples, c(beta.sample))
	r.samples = cbind(r.samples, c(r.sample))
	sigma2.samples= c(sigma2.samples, sigma2.sample)

		iter = iter+1

    if(iter==num.of.iteration)
      if(bm(sigma2.samples)$se > MCSE.Sigma2.Given){
        num.of.iteration = num.of.iteration + 100
        cat("Sigma2 = ", bm(sigma2.samples)$est, "\tMCSE = ", bm(sigma2.samples)$se, "\tNumber of Iterations = ", num.of.iteration, "\n")
      }
	
  }
  end.run = proc.time()
  list(beta.samples = beta.samples, sigma2.samples = sigma2.samples, r.samples = r.samples, Iteration = num.of.iteration, TimeElapsed = (end.run-start.run))
}
CompWiseGibbsSimple = cmpfun(CompWiseGibbsSimple)
