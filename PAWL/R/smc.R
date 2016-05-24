###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################

## utils
fastrmvnorm <- function(n, mu, sigma = diag(length(mu))){
  ev <- eigen(sigma, symmetric = TRUE)
  retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mu, "+")
  return(retval)
}

fastlogdmvnorm <- function(x, mu, Sigma){
  distval <- mahalanobis(x, center = mu, cov = Sigma)
  logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  return(logretval)
}

################################################
## SMC Sampler
################################################

## Takes weights on the log scale and returns normalized weights
normalizeweight <- function(log_weights){
  corrected_log_weights <- log_weights - max(log_weights)
  normalized_weights <- exp(corrected_log_weights) / sum(exp(corrected_log_weights))
  return(normalized_weights)
}
## Multinomial resampling
## Takes normalized weights and returns indices
multinomial_resampling <- function(normalized_weights){
  return(sample(1:length(normalized_weights), length(normalized_weights), replace = TRUE, prob = normalized_weights))
}
## Residual resampling
## Takes normalized weights and returns indices
residual_resampling <- function(normalized_weights){
	l <- length(normalized_weights)
	n_j <- l * normalized_weights
	r_j <- n_j - floor(n_j)
	n_j <- floor(n_j)
	indices <- c()
	for (index in 1:l){
        indices <- c(indices, rep(index, n_j[index]))
	}
	if (sum(r_j) > 0)
        indices <- c(indices, sample(1:l, size= (l - length(indices)), replace = T, prob = r_j))
	return(indices)
}
## Systematic resampling
## Takes normalized weights and returns indices
systematic_resampling <- function(normalized_weights){
  N <- length(normalized_weights)
  indices <- rep(0, N)
  normalized_weights <- N * normalized_weights
  j <- 1
  csw <- normalized_weights[1]
  u <- runif(1, min = 0, max = 1)
  for (k in 1:N){
    while (csw < u){
      j <- j + 1
      csw <- csw + normalized_weights[j]
    }
    indices[k] <- j
    u <- u + 1
  }
  return(indices)
}
## Effective Sample Size
## Takes normalized weights and returns ESS
ESSfunction <- function(normalized_weights){
  sqweights <- normalized_weights**2
  return(1. / sum(sqweights))
}

#### Sequential Monte Carlo sampler on a sequence of tempered target distributions
smc <- function(target, AP, verbose = TRUE){
  if (verbose) print("Launching SMC with parameters:")
  if (verbose) print(AP)
  size <- AP@nparticles
  temperatures <- AP@temperatures
  particles <- as.matrix(target@rinit(size), ncol = target@dimension)
  nbsteps <- length(temperatures)
  nmoves <- AP@nmoves
  targetlogdensity <- function(x) target@logdensity(x, target@parameters)
  logtargetvalues <- targetlogdensity(particles)
  loginitvalues <- target@dinit(particles, target@parameters)
  weights <- (- temperatures[1]) * loginitvalues + temperatures[1] * logtargetvalues
  ESSarray <- rep(0, nbsteps - 1)
  acceptratios <- c()
  resamplingtimes <- c()
  if (AP@movetype == "independent"){
    rproposal <- function(x, muhat, sigmahat){
      return(fastrmvnorm(size, mu = muhat, sigma = sigmahat))
    }
    dproposal <- function(x, muhat, sigmahat) fastlogdmvnorm(x, muhat, sigmahat)
  } else { #random walk
    rproposal <- function(x, muhat, sigmahat){
      return(x + fastrmvnorm(size, mu = rep(0, target@dimension), 
                                     sigma = AP@movescale * sigmahat))
    }
    dproposal <- function(x, muhat, sigmahat) return(0)
  }
  if (AP@resamplingscheme == "multinomial"){
      resampling <- multinomial_resampling
  }
  if (AP@resamplingscheme == "residual"){
      resampling <- residual_resampling
  }
  if (AP@resamplingscheme == "systematic"){
      resampling <- systematic_resampling 
  }
  for (i in 2:nbsteps){
    weights <- weights + (temperatures[i] - temperatures[i-1]) * logtargetvalues +
                         (temperatures[i-1] - temperatures[i]) * loginitvalues
    normweights <- normalizeweight(weights)
    currentESS <- ESSfunction(normweights)
    ESSarray[i-1] <- currentESS
    if (currentESS < AP@ESSthreshold * size){
      if (verbose) cat("Step", i, " reaching temperature ", temperatures[i], "\n")
      if (verbose) cat("ESS (in %):", currentESS / size, "\n")
      if (verbose) cat("***** resample-move\n")
      resamplingtimes <- c(resamplingtimes, i)
      resampled_indices <- resampling(normweights)
      particles <- as.matrix(particles[resampled_indices,], ncol = target@dimension)
      logtargetvalues <- logtargetvalues[resampled_indices]
      loginitvalues <- loginitvalues[resampled_indices]
      weights <- rep(0, size)
      if (nmoves > 0){
        for (moveindex in 1:nmoves){
          Muhat <- apply(particles, 2, mean)
          Sigmahat <- cov(particles)
          if (all(Sigmahat == 0)){
            Sigmahat <- Sigmahat + diag(0.01, rep(target@dimension))
          }
          proposals <- rproposal(particles, Muhat, Sigmahat)
          omegacurrent <- temperatures[i] * logtargetvalues + (1 - temperatures[i]) * loginitvalues - 
                            dproposal(particles, Muhat, Sigmahat)
          targetprop <- targetlogdensity(proposals)
          initprop <- target@dinit(proposals, target@parameters)
          omegaproposals <- temperatures[i] * targetprop + (1 - temperatures[i]) * initprop -
                            dproposal(proposals, Muhat, Sigmahat)
          acceptations <- (log(runif(size)) < (omegaproposals - omegacurrent))
          particles[acceptations,] <- proposals[acceptations,]
          logtargetvalues[acceptations] <- targetprop[acceptations]
          loginitvalues[acceptations] <- initprop[acceptations]
          acceptratio <- mean(acceptations)
          if (verbose) cat("Acceptance rate:", acceptratio, "\n")
        }
      }
    }
  }
  return(list(particles = particles, 
              weights = weights,
              ESSarray = ESSarray,
              resamplingtimes = resamplingtimes))
}

