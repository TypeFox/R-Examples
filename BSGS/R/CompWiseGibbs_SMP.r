#=======================================================================================================================================#
# Bayesian Variable Selection Approach with Gibbs Sampler                                                                               #
# Date: 02/28/2013                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Perform Bayesian variable selection with component-wise prior to identify the important variable                         #
#=======================================================================================================================================#

CompWiseGibbsSMP = function(Y, X, beta.value, r, tau2, rho, sigma2, nu, lambda, num.of.inner.iter, num.of.iteration, MCSE.Sigma2.Given)
{

	num.of.obs = length(Y)
	num.of.covariates = ncol(X)
	x = X
	y = Y
	r.samples = sigma2.samples = beta.samples = numeric(0)

  	r.sample = r.samples = r
  	beta.sample = beta.samples = c(beta.value)
  	sigma2.sample = sigma2.samples = sigma2

	iter = 1
	start.run = proc.time()
	while(iter < num.of.iteration){
	#===========================================
	# Add a variable
	#===========================================
		for(num.of.inner.iter.index in 1:num.of.inner.iter){
			Add.or.Delete = rbinom(1,1, 0.5)

			if(all(r.sample ==1)){
				Add.or.Delete = 0
			}
			if(any(r.sample ==0)){

				SSX = apply((as.matrix(x[, r.sample==0]))^2, 2, sum)

				SSR.fun = function(zero.index) t(y -x[, -zero.index] %*% matrix(beta.sample[-zero.index], ncol=1)) %*% as.matrix( x[, zero.index])

				SSR = mapply(SSR.fun, which(r.sample==0),  SIMPLIFY = TRUE)


				mu.beta = SSR * tau2 / (SSX * tau2 + sigma2.sample)
				var.beta = tau2 * sigma2.sample /(SSX * tau2 + sigma2.sample)

				#cat("mu.beta = ", mu.beta, "\n")
				#cat("var.beta = ", var.beta, "\n")

				log.ratio.1.0 = 0.5 * ( log(var.beta) - log(tau2) ) + 0.5 * mu.beta^2/var.beta

				sum.prob.z.0 = sum(exp(log.ratio.1.0))
				
				if(all(r.sample==0))
					Add.or.Delete = 1
			}

			if(Add.or.Delete){

					num.of.active = sum(r.sample)

					prob.acc.add  = sum.prob.z.0 /(sum(r.sample)+1)

					if(runif(1) < prob.acc.add){
						if(length(which(r.sample==0))==1)
							var.selection = which(r.sample==0)
						else
							var.selection = sample(which(r.sample==0), 1, prob = exp(log.ratio.1.0))
						r.sample[var.selection] = 1
						SSX = sum((x[, var.selection])^2)
						SSR = t(y -x[, -var.selection] %*% matrix(c(beta.sample)[-var.selection], ncol=1)) %*% as.matrix( x[, var.selection])
						mu.beta = SSR * tau2 / (SSX * tau2 + sigma2.sample)
						var.beta = tau2 * sigma2.sample /(SSX * tau2 + sigma2.sample)
						beta.sample[var.selection] = rnorm(1, mu.beta, sqrt(var.beta))
					}
			}

			#===========================================
			# Delete a variable
			#===========================================
			else{
					if(sum(r.sample) == 1)
						var.delete = which(r.sample==1)
					var.delete = var.selection = sample(which(r.sample==1), 1)
					SSX = sum((x[, var.delete])^2)
					SSR = t(y -x[, -var.delete] %*% matrix(c(beta.sample)[-var.delete], ncol=1)) %*% as.matrix( x[, var.delete])

					mu.beta = SSR * tau2 / (SSX * tau2 + sigma2.sample)
					var.beta = tau2 * sigma2.sample /(SSX * tau2 + sigma2.sample)


					if(all(r.sample ==1)) sum.prob.z.0=0
					log.ratio.1.0.delete =0.5 * ( log(var.beta) - log(tau2) ) + 0.5 * mu.beta^2/var.beta

					prob.acc.delete  = sum(r.sample) / (sum.prob.z.0 + exp(log.ratio.1.0.delete) )
					if(runif(1) < prob.acc.delete){
						r.sample[var.delete] = 0
						beta.sample[var.delete] =0
					}
					else
						beta.sample[var.delete] = rnorm(1, mu.beta, sqrt(var.beta))
			}
		}
		
		beta.samples = cbind(beta.samples, c(beta.sample))
		r.samples = cbind(r.samples, c(r.sample))
		
		a.sample = num.of.obs+ sum(r.sample)
		b.sample =  sum( (y -x %*% matrix(beta.sample, ncol=1))^2 ) + nu*lambda

		sigma2.sample = rigamma(1, a.sample/2, b.sample/2)

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

