require(coda)


################################################################################
#
# @brief General Markov chain Monte Carlo algorithm using Metropolis-Hastings.
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-11-06, version 1.1
#
# @param    likelihoodFunction     function      the log-likelihood function
# @param    priors                 list          a list of functions of the log-prior density per parameter
# @param    parameters             vector        initial set of parameter values
# @param    logTransform           vector        should be a log-transform be used for parameter proposals
# @param    iterations             scalar        the number of iterations
# @param    burnin                 scalar        number of burnin iterations
# @param    thinning               scalar        number of iterations between samples
# @param    adaptive               boolean       should we use adaptive MCMC?
# @param    verbose                boolean       should we print information during the MCMC?
# @return                          list          samples from the posterior
#
################################################################################

tess.mcmc <- function(likelihoodFunction,priors,parameters,logTransforms,delta,iterations,burnin=round(iterations/3),thinning=1,adaptive=TRUE,verbose=FALSE) {

  OPTIMIZATIONS <- max(ceiling(burnin/20),50)

  # create a list for the samples
  chain = array(dim = c(floor(iterations/thinning)+1,length(parameters))) #reserve memory for the chain, for large chains we might consider writing to a file instead of storing in memory

  # pre-compute current posterior probability
  pp <- likelihoodFunction(parameters)
  for ( j in 1:length(parameters) ) {
    pp <- pp + priors[[j]](parameters[j])
  }

  if ( verbose ) {
    cat("Burning-in the chain ...\n")
    cat("0--------25--------50--------75--------100\n")
    bar <- txtProgressBar(style=1,width=42)
  }

  # this is the tuning parameter
#  delta <- rep(1,length(parameters))
  accepted <- rep(0,length(parameters))
  tried <- rep(0,length(parameters))
  
  # We need an initial sample if the burnin == 0
  if (burnin == 0) {
     chain[1,] <- parameters    
  }

  for (i in 1:(burnin+iterations)) {

    if ( verbose == TRUE ) {
      if ( i <= burnin ) {
        setTxtProgressBar(bar,i/burnin)
      } else if (i == (burnin+1) ) {
        cat("\nFinished burnin period!\n\n")
        cat("Running the chain ...\n")
        cat("0--------25--------50--------75--------100\n")
        bar <- txtProgressBar(style=1,width=42)
      } else {
        setTxtProgressBar(bar,(i-burnin)/iterations)
      }
    }

    # if we use adaptive MCMC we might have to modify our parameters
    if ( adaptive == TRUE ) {
      if ( i <= burnin ) {
        if ( i %% OPTIMIZATIONS == 0 ) {
          for ( j in 1:length(parameters) ) {
            rate <- accepted[j] / tried[j]
            if ( rate > 0.44 ) {
              delta[j] <- delta[j] * (1.0 + ((rate-0.44)/0.56) )
            } else {
              delta[j] <- delta[j] / (2.0 - rate/0.44 )
            }
            tried[j] <- 0
            accepted[j] <- 0
          }
        }
      }
    }

    # propose new values
    for ( j in 1:length(parameters) ) {
      # increase our tried counter
      tried[j] <- tried[j] + 1

      if ( logTransforms[j] == TRUE ) {
        if (parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }
        eta           <- log(parameters[j]) ### propose a new value for parameter[j]
        new_eta       <- eta + rnorm(1,0,delta[j])
        new_val       <- exp(new_eta)
        hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- 0.0
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        if ( is.finite(new_pp) ) {
          new_pp        <- new_pp + likelihoodFunction(parameters)
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
          accepted[j] <- accepted[j] + 1
        } else {
          parameters[j] <- exp(eta)
        }
      } else {
        eta           <- parameters[j] ### propose a new value for parameter[j]
        new_val       <- eta + rnorm(1,0,delta[j])
        hr            <- 0.0 # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- 0.0
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        if ( is.finite(new_pp) ) {
          new_pp        <- new_pp + likelihoodFunction(parameters)
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
          accepted[j] <- accepted[j] + 1
        } else {
          parameters[j] <- eta
        }

      }

    }

    # sample the parameter
    if (i >= burnin) {
      if ( (i-burnin) %% thinning == 0 ) {
        chain[(i-burnin)/thinning+1,] <- parameters
      }
    }

  }

  if ( is.null(names(priors)) ) {
    names(priors) <- 1:length(priors)
  }

  chain <- as.mcmc(chain)
  varnames(chain) <- names(priors)

  if ( verbose ) {
    cat("\nFinished MCMC!\n\n")
    cat("Parameter | delta | Acceptance Probability\n")
    cat("==========================================\n")
    for ( j in 1:length(parameters) ) {
      cat(sprintf("%s\t\t| %.3f\t\t| %.3f\n",names(priors)[j],delta[j],accepted[j]/tried[j]))
    }
  }

  return(chain) #return a mcmc object, used by coda to plot

}






