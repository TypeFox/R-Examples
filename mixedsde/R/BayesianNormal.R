#' Bayesian Estimation In Mixed Stochastic Differential Equations
#' 
#' @description Gibbs sampler for Bayesian estimation of the random effects \eqn{(\alpha_j, \beta_j)} in the mixed SDE
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma a(X_j(t)) dW_j(t)}.
#' @param times vector of observation times
#' @param X matrix of the M trajectories (each row is a trajectory with \eqn{N= T/\Delta} column).
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param prior list of prior parameters: mean and variance of the Gaussian prior on the mean mu, shape and scale of the inverse Gamma prior for the variances omega, shape and scale of the inverse Gamma prior for sigma  
#' @param start list of starting values: mu, sigma  
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param nMCMC number of iterations of the MCMC algorithm
#' @param propSd proposal standard deviation of \eqn{\phi} is \eqn{|\mu|*}\code{propSd}\eqn{/\log(N)} at the beginning, is adjusted when acceptance rate is under 30\% or over 60\%
#' @return
#' \item{alpha}{posterior samples (Markov chain) of \eqn{\alpha}}
#' \item{beta}{posterior samples (Markov chain) of \eqn{\beta}}
#' \item{mu}{posterior samples (Markov chain) of \eqn{\mu}}
#' \item{omega}{posterior samples (Markov chain) of \eqn{\Omega}}
#' \item{sigma2}{posterior samples (Markov chain) of \eqn{\sigma^2}}
#' @references 
#' Hermann, S., Ickstadt, K. and C. Mueller (2016). Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. \emph{Appearing in: Applied Stochastic Models in Business and Industry}. 
#' 
#' Rosenthal, J. S. (2011). 'Optimal proposal distributions and adaptive MCMC.' Handbook of Markov Chain Monte Carlo (2011): 93-112.

BayesianNormal <- function(times, X, model = c("OU", "CIR"), prior, start, random, nMCMC = 1000, propSd = 0.2) {
    model <- match.arg(model)
    M <- nrow(X)
    N <- ncol(X)
    delta <- diff(times)
    propSdSigma2 <- 0.01
    propSdPhi <- c(abs(start$mu) * propSd/log(N))
    
    if (model == "OU") {
        likeli <- dcOU
        postSigma <- function(phi, sigmaOld = 1, propSdSigma2 = NULL) {
            
            alphaPost <- prior$alpha.sigma + M * (N - 1)/2
            
            help <- numeric(M)
            for (i in 1:M) {
                help[i] <- sum(phi[i, 2]/(1 - exp(-2 * phi[i, 2] * delta)) * (X[i, -1] - phi[i, 1]/phi[i, 2] - (X[i, -N] - phi[i, 1]/phi[i, 
                  2]) * exp(-phi[i, 2] * delta))^2)
            }
            betaPost <- prior$beta.sigma + sum(help)
            1/rgamma(1, alphaPost, betaPost)
        }
    } else {
        likeli <- dcCIR2
        postSigma <- function(phi, sigmaOld = 1, propSdSigma2 = NULL) {
            # alphaPost <- prior$alpha.sigma + M*(N-1)/2 help <- numeric(M) for (i in 1:M) { help[i] <- sum((X[i,-1] - X[i,-N] - b(phi[i, ],
            # times[-N], X[i,-N]) * delta)^2/ (X[i,-N] * delta)) } betaPost <- prior$beta.sigma + sum(help)/2 return(1/rgamma(1, alphaPost,
            # betaPost))
            
            sigma_drawn <- rnorm(1, sqrt(sigmaOld), propSdSigma2)^2
            ratio <- dgamma(1/sigma_drawn, prior$alpha.sigma, prior$beta.sigma)/dgamma(1/sigmaOld, prior$alpha.sigma, prior$beta.sigma)  # prior ratio
            ratio <- ratio * sqrt(sigma_drawn)/sqrt(sigmaOld)  # proposal ratio
            he <- sapply(1:M, function(i) {
                prod(dcCIR2(X[i, -1], delta, X[i, -N], c(phi[i, ], sqrt(sigma_drawn)))/dcCIR2(X[i, -1], delta, X[i, -N], c(phi[i, ], sqrt(sigmaOld))))
            })
            ratio <- ratio * prod(he)
            if (is.na(ratio)) 
                ratio <- 0
            if (runif(1) <= ratio) {
                return(sigma_drawn)
            } else {
                return(sigmaOld)
            }
        }
    }
    if (length(random) == 2) {
        
        postPhij <- function(lastPhi, mu, Omega, sigma, Xj, propSdPhi) {
            phi_out <- lastPhi
            phi <- lastPhi
            for (k in 1:2) {
                phi[k] <- phi_out[k] + rnorm(1, 0, propSdPhi[k])
                ratio <- prod(dnorm(phi, mu, sqrt(Omega))/dnorm(phi_out, mu, sqrt(Omega)))
                
                ratio <- ratio * prod(likeli(Xj[-1], delta, Xj[-N], c(phi, sqrt(sigma)))/likeli(Xj[-1], delta, Xj[-N], c(phi_out, sqrt(sigma))))
                # ratio <- ratio * prod(dnorm(Xj[-1], Xj[-N] + b(phi, times[-N], Xj[-N]) * delta, sqrt(sigma * Xj[-N] * delta))/ dnorm(Xj[-1],Xj[-N] +
                # b(phi_out, times[-N], Xj[-N]) * delta, sqrt(sigma * Xj[-N] * delta)))
                if (is.na(ratio)) {
                  ratio <- 0
                }
                if (runif(1) < ratio) {
                  phi_out[k] <- phi[k]
                } else {
                }
            }
            phi_out
        }
        postmu <- function(phi, Omega) {
            Vpost <- 1/(1/prior$v + M/Omega)
            mpost <- Vpost * ((1/prior$v) %*% prior$m + apply((1/Omega) * t(phi), 1, sum))
            
            rnorm(2, mpost, sqrt(Vpost))
        }
        postOm <- function(phi, mu) {
            Dia <- numeric(2)
            for (i in 1:2) {
                Dia[i] <- 1/rgamma(1, prior$alpha.omega[i] + M/2, prior$beta.omega[i] + sum((phi[, i] - mu[i])^2)/2)
            }
            Dia
        }
        # storage and starting variables
        alpha_out <- matrix(0, nMCMC, M)
        beta_out <- matrix(0, nMCMC, M)
        mu_out <- matrix(0, nMCMC, length(start$mu))
        Omega_out <- matrix(0, nMCMC, length(start$mu))
        phi <- matrix(rep(start$mu, each = M), M)
        mu <- start$mu
        Omega <- postOm(phi, mu) + 0.1
        
    } else {
        postRandom <- function(lastPhi, mu, Omega, sigma, Xj, propSdRandom) {
            phi_out <- lastPhi
            phi <- lastPhi
            phi[random] <- phi[random] + rnorm(1, 0, propSdRandom)
            ratio <- dnorm(phi[random], mu, sqrt(Omega))/dnorm(phi_out[random], mu, sqrt(Omega))
            ratio <- ratio * prod(likeli(Xj[-1], delta, Xj[-N], c(phi, sqrt(sigma)))/likeli(Xj[-1], delta, Xj[-N], c(phi_out, sqrt(sigma))))
            # ratio <- ratio * prod(dnorm(Xj[-1], Xj[-N] + b(phi, times[-N], Xj[-N]) * delta, sqrt(sigma * Xj[-N] * delta))/ dnorm(Xj[-1], Xj[-N]
            # + b(phi_out, times[-N], Xj[-N]) * delta, sqrt(sigma * Xj[-N] * delta)))
            if (is.na(ratio)) {
                ratio <- 0
            }
            if (runif(1) < ratio) {
                phi_out[random] <- phi[random]
            } else {
            }
            phi_out[random]
        }
        postmu <- function(phi, Omega) {
            n <- length(phi)
            Vpost <- 1/(1/prior$v[random] + n * 1/Omega)
            mpost <- Vpost %*% ((1/prior$v[random]) * prior$m[random] + sum(1/Omega * phi))
            
            rnorm(1, mpost, sqrt(Vpost))
        }
        postOm <- function(phi, mu) {
            1/rgamma(1, prior$alpha.omega + length(phi)/2, prior$beta.omega + sum((phi - mu)^2)/2)
        }
        # storage and starting variables
        mu_out <- numeric(nMCMC)
        Omega_out <- numeric(nMCMC)
        
        if (random == 1) {
            postNonrandom <- function(lastRandom, lastNonrandom, sigma, propSdFixed) {
                
                beta <- lastNonrandom
                beta_drawn <- beta + rnorm(1, 0, propSdFixed)
                ratio <- dnorm(beta_drawn, prior$m[2], prior$v[2])/dnorm(beta, prior$m[2], prior$v[2])
                he <- sapply(1:M, function(i) {
                  prod(likeli(X[i, -1], delta, X[i, -N], c(lastRandom[i], beta_drawn, sqrt(sigma)))/likeli(X[i, -1], delta, X[i, -N], 
                    c(lastRandom[i], beta, sqrt(sigma))))
                })
                ratio <- ratio * prod(he)
                if (is.na(ratio)) {
                  ratio <- 0
                }
                if (runif(1) <= ratio) {
                  beta <- beta_drawn
                } else {
                }
                beta
            }
            # storage and starting variables
            alpha_out <- matrix(0, nMCMC, M)
            beta_out <- numeric(nMCMC)
            alpha <- rep(start$mu[1], M)
            beta <- start$mu[2]
            mu <- start$mu[1]
            Omega <- postOm(alpha, mu) + 0.1
        }
        if (random == 2) {
            postNonrandom <- function(lastRandom, lastNonrandom, sigma, propSdFixed) {
                
                alpha <- lastNonrandom
                alpha_drawn <- alpha + rnorm(1, 0, propSdFixed)
                ratio <- dnorm(alpha_drawn, prior$m[1], prior$v[1])/dnorm(alpha, prior$m[1], prior$v[1])
                he <- sapply(1:M, function(i) {
                  prod(likeli(X[i, -1], delta, X[i, -N], c(alpha_drawn, lastRandom[i], sqrt(sigma)))/likeli(X[i, -1], delta, X[i, -N], 
                    c(alpha, lastRandom[i], sqrt(sigma))))
                })
                # he <- sapply(1:M, function(i){ prod(dnorm(X[i,-1], X[i,-N] + b(c(alpha_drawn, lastRandom[i]), times[-N], X[i,-N]) * delta,
                # sqrt(sigma * X[i,-N] * delta))/ dnorm(X[i,-1], X[i,-N] + b(c(alpha, lastRandom[i]), times[-N], X[i,-N]) * delta, sqrt(sigma *
                # X[i,-N] * delta))) } )
                ratio <- ratio * prod(he)
                if (is.na(ratio)) {
                  ratio <- 0
                }
                if (runif(1) <= ratio) {
                  alpha <- alpha_drawn
                } else {
                }
                alpha
            }
            # storage and starting variables
            alpha_out <- numeric(nMCMC)
            beta_out <- matrix(0, nMCMC, M)
            alpha <- start$mu[1]
            beta <- rep(start$mu[2], M)
            mu <- start$mu[2]
            Omega <- postOm(beta, mu) + 0.1
        }
    }
    # storage and starting variables
    sigma_out <- numeric(nMCMC)
    sigma <- start$sigma
    
    if (length(random) == 2) {
        
        for (count in 1:nMCMC) {
            for (i in 1:M) {
                phi[i, ] <- postPhij(phi[i, ], mu, Omega, sigma, X[i, ], propSdPhi)
            }
            
            mu <- postmu(phi, Omega)
            Omega <- postOm(phi, mu)
            
            sigma <- postSigma(phi, sigma, propSdSigma2)
            
            alpha_out[count, ] <- phi[, 1]
            beta_out[count, ] <- phi[, 2]
            mu_out[count, ] <- mu
            Omega_out[count, ] <- Omega
            sigma_out[count] <- sigma
            if (count%%1000 == 0) {
                message(paste(count, "iterations done"))
            }
            
            if (count%%50 == 0) {
                propSdPhi <- c(ad.propSd_random(alpha_out[(count - 50 + 1):count, ], propSdPhi[1], count), ad.propSd_random(beta_out[(count - 
                  50 + 1):count, ], propSdPhi[2], count/50))
            }
            if (model == "CIR" && count%%50 == 0) {
                propSdSigma2 <- ad.propSd(sigma_out[(count - 50 + 1):count], propSdSigma2, count/50)
            }
        }
    } else {
        if (random == 1) {
            for (count in 1:nMCMC) {
                for (i in 1:M) {
                  alpha[i] <- postRandom(c(alpha[i], beta), mu, Omega, sigma, X[i, ], propSdPhi[random])
                }
                beta <- postNonrandom(alpha, beta, sigma, propSdPhi[-random])
                
                mu <- postmu(alpha, Omega)
                Omega <- postOm(alpha, mu)
                
                sigma <- postSigma(cbind(alpha, beta), sigma, propSdSigma2)
                
                alpha_out[count, ] <- alpha
                beta_out[count] <- beta
                mu_out[count] <- mu
                Omega_out[count] <- Omega
                sigma_out[count] <- sigma
                if (count%%1000 == 0) {
                  message(paste(count, "iterations done"))
                }
                
                
                if (count%%50 == 0) {
                  propSdPhi <- c(ad.propSd_random(alpha_out[(count - 50 + 1):count, ], propSdPhi[1], count), ad.propSd(beta_out[(count - 
                    50 + 1):count], propSdPhi[2], count/50))
                }
                
                if (model == "CIR" && count%%50 == 0) {
                  propSdSigma2 <- ad.propSd(sigma_out[(count - 50 + 1):count], propSdSigma2, count/50)
                }
                
            }
        }
        if (random == 2) {
            
            for (count in 1:nMCMC) {
                for (i in 1:M) {
                  beta[i] <- postRandom(c(alpha, beta[i]), mu, Omega, sigma, X[i, ], propSdPhi[random])
                }
                alpha <- postNonrandom(beta, alpha, sigma, propSdPhi[-random])
                
                mu <- postmu(beta, Omega)
                Omega <- postOm(beta, mu)
                
                sigma <- postSigma(cbind(alpha, beta), sigma, propSdSigma2)
                
                alpha_out[count] <- alpha
                beta_out[count, ] <- beta
                mu_out[count] <- mu
                Omega_out[count] <- Omega
                sigma_out[count] <- sigma
                if (count%%1000 == 0) {
                  message(paste(count, "iterations done"))
                }
                
                if (count%%50 == 0) {
                  propSdPhi <- c(ad.propSd(alpha_out[(count - 50 + 1):count], propSdPhi[1], count), ad.propSd_random(beta_out[(count - 
                    50 + 1):count, ], propSdPhi[2], count/50))
                }
                
                if (model == "CIR" && count%%50 == 0) {
                  propSdSigma2 <- ad.propSd(sigma_out[(count - 50 + 1):count], propSdSigma2, count/50)
                }
                
            }
        }
    }
    
    list(alpha = alpha_out, beta = beta_out, mu = mu_out, omega = Omega_out, sigma2 = sigma_out)
}

######## 
#' Likelihood Function For The CIR Model
#' 
#' @description Likelihood
#' @param x current observation
#' @param t time of observation
#' @param x0 starting point, i.e. observation in time 0
#' @param theta parameter \eqn{(\alpha, \beta, \sigma)}
#' @param log logical(1) if TRUE, log likelihood
#' @references 
#' Iacus, S. M. (2008). Simulation and Inference for Stochastic Differential Equations.
#' 
dcCIR2 <- function(x, t, x0, theta, log = FALSE) {
    c <- 2 * theta[2]/((1 - exp(-theta[2] * t)) * theta[3]^2)
    ncp <- 2 * c * x0 * exp(-theta[2] * t)
    df <- 4 * theta[1]/theta[3]^2
    lik <- (dchisq(2 * x * c, df = df, ncp = ncp, log = TRUE) + log(2 * c))
    if (!log) 
        lik <- exp(lik)
    lik
}
# the dcCIR function from the sde package does not work...

######## 
#' Removing Of Burn-in Phase And Thinning
#' 
#' @description Transfers class object Bayes.fit from the original to the thinned chains
#' @param res Bayes.fit class object
#' @param burnIn number of burn-in samples
#' @param thinning thinning rate
chain2samples <- function(res, burnIn, thinning) {
    if (missing(burnIn)) 
        burnIn <- res@burnIn
    if (missing(thinning)) 
        thinning <- res@thinning
    
    ind.chain <- seq(burnIn + 1, length(res@sigma2), by = thinning)
    return(new(Class = class(res), prior = res@prior, alpha = as.matrix(res@alpha[ind.chain, ]), beta = as.matrix(res@beta[ind.chain, 
        ]), random = res@random, mu = as.matrix(res@mu[ind.chain, ]), omega = as.matrix(res@omega[ind.chain, ]), sigma2 = res@sigma2[ind.chain], 
        model = res@model, times = res@times, X = res@X))
    
}

######## 
#' Calcucation Of Burn-in Phase And Thinning Rate
#' 
#' @description Proposal for burn-in and thin rate
#' @param results Bayes.fit class object
#' @param random one out of 1, 2, c(1,2)
diagnostic <- function(results, random) {
    diagnostic.own <- function(chain, dependence = 0.8, m = 10) {
        lc <- length(chain)
        K <- floor(lc/m)
        thinning <- min(which(acf(chain[-(1:floor(lc/5))], plot = F)$acf <= dependence)[1], floor(K/10), na.rm = TRUE)
        he1 <- sapply(1:m, function(i) chain[((i - 1) * K + 1):(i * K)][seq(1, K, by = thinning)])
        
        he2 <- apply(he1, 2, quantile, c(0.025, 0.975))
        he.mean <- apply(he1, 2, mean)
        is.in <- (he.mean[-m] >= he2[1, -1] & he.mean[-m] <= he2[2, -1]) | (he.mean[-1] >= he2[1, -m] & he.mean[-1] <= he2[2, -m])
        # burnIn <- 0
        burnIn <- K
        for (i in 1:(m - 1)) {
            if (sum(is.in) < length(is.in)) {
                is.in <- is.in[-1]
                burnIn <- burnIn + K
            }
        }
        burnIn <- min((m - 1) * K, burnIn)
        return(c(burnIn = burnIn, thinning = thinning))
    }
    if (length(random) == 2) {
        # he <- matrix(0,5,4) he[1, ] <- raftery.diag(as.mcmc(results$mu[,1]))$resmatrix he[2, ] <-
        # raftery.diag(as.mcmc(results$mu[,2]))$resmatrix he[3, ] <- raftery.diag(as.mcmc(results$omega[,1]))$resmatrix he[4, ] <-
        # raftery.diag(as.mcmc(results$omega[,2]))$resmatrix he[5, ] <- raftery.diag(as.mcmc(results$sigma2))$resmatrix burnIn <- max(he[,1])
        # thinning <- floor(length(results$sigma2)/max(he[,3]))
        he <- matrix(0, 5, 2)
        he[1, ] <- diagnostic.own(results$mu[, 1])
        he[2, ] <- diagnostic.own(results$mu[, 2])
        he[3, ] <- diagnostic.own(results$omega[, 1])
        he[4, ] <- diagnostic.own(results$omega[, 2])
        he[5, ] <- diagnostic.own(results$sigma2)
        burnIn <- max(he[, 1])
        thinning <- max(he[, 2])
        
    } else {
        if (random == 1) {
            # he <- matrix(0,4,4) he[1, ] <- raftery.diag(as.mcmc(results$mu))$resmatrix he[2, ] <- raftery.diag(as.mcmc(results$beta))$resmatrix
            # he[3, ] <- raftery.diag(as.mcmc(results$omega))$resmatrix he[4, ] <- raftery.diag(as.mcmc(results$sigma2))$resmatrix burnIn <-
            # max(he[,1]) thinning <- floor(length(results$sigma2)/max(he[,3]))
            he <- matrix(0, 4, 2)
            he[1, ] <- diagnostic.own(results$mu)
            he[2, ] <- diagnostic.own(results$beta)
            he[3, ] <- diagnostic.own(results$omega)
            he[4, ] <- diagnostic.own(results$sigma2)
            burnIn <- max(he[, 1])
            thinning <- max(he[, 2])
            
        } else {
            # he <- matrix(0,4,4) he[1, ] <- raftery.diag(as.mcmc(results$mu))$resmatrix he[2, ] <- raftery.diag(as.mcmc(results$alpha))$resmatrix
            # he[3, ] <- raftery.diag(as.mcmc(results$omega))$resmatrix he[4, ] <- raftery.diag(as.mcmc(results$sigma2))$resmatrix burnIn <-
            # max(he[,1]) thinning <- floor(length(results$sigma2)/max(he[,3]))
            he <- matrix(0, 4, 2)
            he[1, ] <- diagnostic.own(results$mu)
            he[2, ] <- diagnostic.own(results$alpha)
            he[3, ] <- diagnostic.own(results$omega)
            he[4, ] <- diagnostic.own(results$sigma2)
            burnIn <- max(he[, 1])
            thinning <- max(he[, 2])
            
        }
    }
    
    thinning <- min(thinning, ceiling((length(results$sigma2) - burnIn)/100))
    
    return(list(burnIn = burnIn, thinning = thinning))
}



#' Adaptation For The Proposal Variance
#' 
#' @description Calculation of new proposal standard deviation
#' @param chain vector of Markov chain samples
#' @param propSd old proposal standard deviation
#' @param iteration number of current iteration
#' @param lower lower bound
#' @param upper upper bound
#' @param delta.n function for adding/subtracting from the log propSd
#' @references 
#' Rosenthal, J. S. (2011). Optimal proposal distributions and adaptive MCMC. Handbook of Markov Chain Monte Carlo, 93-112.
ad.propSd <- function(chain, propSd, iteration, lower = 0.3, upper = 0.6, delta.n = function(n) min(0.1, 1/sqrt(n))) {
    ar <- length(unique(chain))/length(chain)
    lsi <- log(propSd)
    
    lsi[ar < lower] <- lsi - delta.n(iteration)
    lsi[ar > upper] <- lsi + delta.n(iteration)
    exp(lsi)
    
}

#' Adaptation For The Proposal Variance
#' 
#' @description Calculation of new proposal standard deviation for the random effects
#' @param chain matrix of Markov chain samples
#' @param propSd old proposal standard deviation
#' @param iteration number of current iteration
#' @param lower lower bound
#' @param upper upper bound
#' @param delta.n function for adding/subtracting from the log propSd
#' @references 
#' Rosenthal, J. S. (2011). Optimal proposal distributions and adaptive MCMC. Handbook of Markov Chain Monte Carlo, 93-112.
ad.propSd_random <- function(chain, propSd, iteration, lower = 0.3, upper = 0.6, delta.n = function(n) min(0.1, 1/sqrt(n))) {
    ar <- apply(chain, 2, function(vec) length(unique(vec))/length(vec))
    lsi <- log(propSd)
    
    lsi[mean(ar) < lower] <- lsi - delta.n(iteration)
    lsi[mean(ar) > upper] <- lsi + delta.n(iteration)
    exp(lsi)
} 
