require(coda)


################################################################################
# 
# @brief General stepping stone sampling to compute marginal likelihood of a model
#
# @date Last modified: 2012-11-09
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-11-07, version 1.0
# @citation Xie et al., 2011, Improving Marginal Likelihood Estimation for Bayesian Phylogenetic Model Selection
#
# @param    likelihoodFunction     function      the likelihood function of the model
# @param    priorFunction          function      the prior function of the model
# @param    parameters             vector        initial set of parameter values
# @param    logTransform           vector        should be a log-transform be used for parameter proposals
# @param    iterations             scalar        the number of iterations
# @param    burnin                 scalar        the number of iterations that should be burned before starting
# @param    K                      scalar        the number of power posteriors
# @return                          scaler        marginal likelihood of the model
#
################################################################################


tess.steppingStoneSampling <- function(likelihoodFunction,priors,parameters,logTransforms,iterations,burnin=round(iterations/3),K=50) {

  x <- K:0 / K
  beta <- qbeta(x,0.3,1)

  # pre-compute current posterior probability
  pp    <- likelihoodFunction(parameters)
  prior <- 0
  for ( j in 1:length(parameters) ) {
    prior <- prior + priors[[j]](parameters[j])
  }

  # the path values
  BF <- 0

  for (k in 1:length(beta)) {
    b <- beta[k]

    samples <- c()
    for (i in 1:(iterations+burnin)) {
    
      # propose new values
      for ( j in 1:length(parameters) ) {
        new_prior <- prior - priors[[j]](parameters[j])
        if ( logTransforms[j] == TRUE ) {
          if (parameters[j] == 0) {
            stop("Cannot propose new value for a parameter with value 0.0.")
          }
          eta           <- log(parameters[j]) ### propose a new value for parameter[j]
          new_eta       <- eta + rnorm(1,0,1)
          new_val       <- exp(new_eta)
          hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
          parameters[j] <- new_val
          new_pp        <- likelihoodFunction(parameters)
          new_prior     <- new_prior + priors[[j]](parameters[j])
        
          # compute acceptance ratio for power posterior
          if ( b == 0.0 ) {
            acceptance_ratio <- new_prior-prior+hr
          } else { 
            acceptance_ratio <- new_prior-prior+b*new_pp-b*pp+hr
          }
          # accept / reject
          if (is.finite(new_pp) && is.finite(new_prior) && acceptance_ratio > log(runif(1,0,1)) ) {
            pp <- new_pp
            prior <- new_prior
          } else {
            parameters[j] <- exp(eta)
          }
        } else {
          eta           <- parameters[j] ### propose a new value for parameter[j]
          new_val       <- eta + rnorm(1,0,1)
          hr            <- 0.0 # calculate the Hastings ratio
          parameters[j] <- new_val
          new_pp        <- likelihoodFunction(parameters)
          new_prior     <- new_prior + priors[[j]](parameters[j])
          # compute acceptance ratio for power posterior
          if ( b == 0.0 ) {
            acceptance_ratio <- new_prior-prior+hr
          } else { 
            acceptance_ratio <- new_prior-prior+b*new_pp-b*pp+hr
          }
          # accept / reject
          if (is.finite(new_pp) && is.finite(new_prior) && acceptance_ratio > log(runif(1,0,1)) ) {
            pp <- new_pp
            prior <- new_prior
          } else {
            parameters[j] <- eta
          }
        }

      }


      # sample the likelihood
      if (i > burnin) {
        samples[i-burnin] <- pp
      }
    }

    max <- samples[which.max(samples)]
    if ( k > 1 ) {
      BF <- BF + log(mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) ))+(beta[k-1]-beta[k])*max
    }
                              
  }

  return (BF) #return the Bayes factor
}
