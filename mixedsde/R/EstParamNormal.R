#' Maximization Of The Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Maximization of the loglikelihood of the mixed SDE with Normal distribution of the random effects
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}, done with \code{\link[=mixedsde]{likelihoodNormal}}
#' @param U matrix of M sufficient statistics U
#' @param V list of the M sufficient statistics matrix V
#' @param K number of times of observations
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param estim.fix 1 if the fixed parameter is estimated, when random 1 or 2 , 0 otherwise
#' @param fixed value of the fixed parameter if known (not estimated)
#' @return
#' \item{mu}{estimated value of the mean}
#' \item{Omega}{estimated value of the variance}
#' \item{BIChere}{BIC indicator}
#' \item{AIChere}{AIC indicator}

EstParamNormal = function(U, V, K, random, estim.fix, fixed = 0) {
    
    if (length(random) == 1) {
        
        
        if (estim.fix == 0) {
            
            ln = function(param) {
                likelihoodNormal(param[1], param[2], Ur, Vrr, estimphi, random)
            }
            M <- dim(U)[2]
            fix <- (random == 1) + 1
            Ur <- rep(0, M)
            Uf <- rep(0, M)
            Vrr <- rep(0, M)
            Vff <- rep(0, M)
            for (j in 1:M) {
                Ur[j] <- U[random, j] - fixed * V[[j]][1, 2]
                Uf[j] <- U[fix, j]
                Vrr[j] <- V[[j]][random, random]
                Vff[j] <- V[[j]][fix, fix]
            }
            estimphi <- Ur/Vrr
            init.mu <- mean(estimphi)
            init.omega <- sd(estimphi)
            
            res <- optim(c(init.mu, init.omega), fn = ln, method = "Nelder-Mead")
            mu <- res$par[1]
            omega <- abs(res$par[2])
            
            BIChere <- ln(c(mu, omega)) + sum(-2 * fixed * Uf + fixed^2 * Vff) + log(length(estimphi)) * 2
            AIChere <- ln(c(mu, omega)) + sum(-2 * fixed * Uf + fixed^2 * Vff) + 2
        }

        if (estim.fix == 1) {
            
            ln = function(param) {
                likelihoodNormalestimfix(param[1], param[2], param[3], U, V, estimphi, random)  
            }
            estimphi <- matrix(0, length(V), 2)
            for (j in 1:length(V)) {
                estimphi[j, ] <- (1/det(V[[j]])) * matrix(c(V[[j]][2, 2], -V[[j]][1, 2], -V[[j]][1, 2], V[[j]][1, 1]), 2, 2) %*% U[, j]
            }
            estimphi <- t(estimphi)
            init.mu <- c(mean(estimphi[1, ]), mean(estimphi[2, ]))
            init.omega <- sd(estimphi[random, ])
            
            res <- optim(c(init.mu, init.omega), fn = ln, method = "Nelder-Mead")
            
            fix <- (random == 1) + 1
            mu <- c(res$par[random], res$par[fix]) * (random == 1) + c(res$par[fix], res$par[random]) * (random == 2)
            omega <- abs(res$par[3])
            
            BIChere <- ln(c(mu, omega)) + log(length(V)) * 2 + log(length(V) * K)
            AIChere <- ln(c(mu, omega)) + 3
            
            
        }
    }
    if (length(random) == 2) {
        
        loglik_Omegadiag <- function(param) {
            mu <- param[1:2]
            omega <- param[3:4]
            
            return(likelihoodNormal(mu, omega, U, V, estimphi, random))
        }
        
        estimphi <- matrix(0, length(V), 2)
        for (j in 1:length(V)) {
            estimphi[j, ] <- (1/det(V[[j]])) * matrix(c(V[[j]][2, 2], -V[[j]][1, 2], -V[[j]][1, 2], V[[j]][1, 1]), 2, 2) %*% U[, j]
        }
        estimphi <- t(estimphi)
        
        init.mu <- c(mean(estimphi[1, ]), mean(estimphi[2, ]))
        init.omega <- c(sd(estimphi[1, ]), sd(estimphi[2, ]))
        
        res <- optim(c(init.mu, init.omega), loglik_Omegadiag, gr = NULL, method = "Nelder-Mead")
        mu <- c(res$par[1], res$par[2])
        omega <- abs(c(res$par[3], res$par[4]))
        BIChere <- loglik_Omegadiag(c(mu, omega)) + log(dim(estimphi)[2]) * 4
        AIChere <- loglik_Omegadiag(c(mu, omega)) + 4
        
    }
    return(list(mu = mu, omega = omega, BIChere = BIChere, AIChere = AIChere))
}

