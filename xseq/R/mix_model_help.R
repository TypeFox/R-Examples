# Some helper functions for mixture model analysis
# 
# Date: 
#   Revised: February 15, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#

#==============================================================================
MixGaussPriorLoglik = function(lambda, mu, sigma, prior, sigma.equal) {
  # Compute the log-likelihood introduced by the priors
  # 
  # Helper function
  
  if(sigma.equal == TRUE) {
    prior.like = -(prior$dof + 2 + length(mu)) * log(sigma[1]) - 
      ((prior$sigma[1])^2 + sum(prior$kapp*(mu - prior$mu)^2))/2/sigma[1]^2
  } else {
    prior.like = -(prior$dof + 3) * log(sigma) - 
      (1 * (prior$sigma)^2 + prior$kapp*(mu - prior$mu)^2)/2/sigma^2
  }

  prior.like = sum(prior.like) + sum(log(lambda) * (prior$alpha - 1))
  
  return(prior.like)
}

MixStudentPriorLoglik = function(lambda, mu, sigma, prior, sigma.equal) {
  prior.like = MixGaussPriorLoglik(lambda, mu, sigma, prior, sigma.equal)
  return(prior.like)
}

#==============================================================================
MixGaussCheckConverge = function(lambda, mu, sigma, loglik.hist, iter, 
                                 max.iter, conv.tol, lambda.tol, verbose) {
  # To check the convergence of EM iterations
  # 
  # Helper function
  
  done = FALSE
  
  loglik.incr = loglik.hist[iter] - loglik.hist[iter-1]
  if (verbose) {
    cat("iteration = ", iter, " loglik increasing = ", loglik.incr, 
        " loglik = ", loglik.hist[iter], "\n" )
    print(rbind(lambda, mu, sigma))
  }
  
  loglik.incr = loglik.incr / abs(loglik.hist[iter])
  if (abs(loglik.incr) < conv.tol) {
    done = TRUE
    return(done)
  }
  
  # Numeric issues
  if(loglik.incr < 0) {
    cat("WARNING! The likelihood decreases", "\n")
    done = TRUE
  }
  
  if(any(lambda < lambda.tol)) {
    cat(paste("WARNING! Mixture components < ", lambda.tol, sep=""), "\n")
    done = TRUE
  }
  
  if (iter == max.iter) {
    cat("WARNING! NOT CONVERGE!", "\n")
    done = TRUE
  }
  
  return(done)
}

#==============================================================================
MixStudentCheckConverge = function(lambda, mu, sigma, nu, loglik.hist, iter, 
                                 max.iter, conv.tol, lambda.tol, verbose) {
  # To check the convergence of EM iterations
  # 
  # Helper function
  
  done = FALSE
  
  loglik.incr = loglik.hist[iter] - loglik.hist[iter-1]
  if (verbose) {
    cat("iteration = ", iter, " loglik increasing = ", loglik.incr, 
        " loglik = ", loglik.hist[iter], "\n" )
    print(rbind(lambda, mu, sigma, nu))
  }
  
  loglik.incr = loglik.incr / abs(loglik.hist[iter])
  if (abs(loglik.incr) < conv.tol) {
    done = TRUE
    return(done)
  }
  
  if(loglik.incr < 0) {
    cat("WARNING! The likelihood decreases", "\n")
    done = TRUE
  }
  
  if(any(lambda < lambda.tol)) {
    cat(paste("WARNING! Mixture components < ", lambda.tol, sep=""), "\n")
    done = TRUE
  }
  
  if (iter == max.iter) {
    cat("WARNING! NOT CONVERGE!", "\n")
    done = TRUE
  }
  
  return(done)
}



