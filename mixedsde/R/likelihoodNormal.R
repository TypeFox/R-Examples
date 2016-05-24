#' Computation Of The Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Computation of -2 loglikelihood of the mixed SDE with Normal distribution of the random effects
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}.
#' @param mu current value of the mean of the normal distribution
#' @param omega current value of the standard deviation of the normal distribution
#' @param U vector of the M sufficient statistics U (see \code{\link{UV}})
#' @param V vector of the M sufficient statistics V (see \code{\link{UV}})
#' @param estimphi vector or matrix of estimators of the random effects 
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}
#' 

likelihoodNormal <- function(mu, omega, U, V, estimphi, random) {
    
    if (length(random) == 1) {
        Omega <- omega^2
        L <- sum(log(1 + Omega * V)) + sum(V/(1 + Omega * V) * (mu - estimphi)^2 - U * estimphi)
    }
    
    
    if (length(random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
        M <- dim(U)[2]
        
        loglik <- vector(length = M)
        I2 <- diag(c(1, 1))
        for (j in 1:M) {
            A <- (I2 + V[[j]] %*% Omega)
            Rinv <- solve(A) %*% V[[j]]
            b <- mu - estimphi[, j]
            loglik[j] <- log(det(A)) + t(b) %*% Rinv %*% b - t(U[, j]) %*% estimphi[, j]
        }
        L <- sum(loglik)
    }
    return(L)
}





#' Likelihood Function When The Fixed Effect Is Estimated
#' 
#' @description Computation of -2 loglikelihood of the mixed SDE with Normal distribution of the random effects when the fixed effect is estimated for random 1 or 2
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}.
#' @param mu1 current value of the mean of the first effect
#' @param mu2 current value of the mean of the second effect
#' @param omega current value of the standard deviation of the normal distribution
#' @param U vector of the M sufficient statistics U (see \code{\link{UV}})
#' @param V vector of the M sufficient statistics V (see \code{\link{UV}})
#' @param estimphi vector or matrix of estimators of the random effects 
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}
#' 


likelihoodNormalestimfix <- function(mu1, mu2, omega, U, V, estimphi, random) {
    if (length(random) == 1) {
        mu <- c(mu1, mu2)
        
        fix <- (random == 1) + 1
        
        psi <- mu[fix]
        
        muRandom <- mu[random]
        
        Omega <- omega^2
        
        V11 <- unlist(V)[seq(1, 4 * length(V), by = 4)]
        V22 <- unlist(V)[seq(4, 4 * length(V), by = 4)]
        
        Vrr <- V11 * (random == 1) + V22 * (random == 2)
        
        Vff <- V11 * (fix == 1) + V22 * (fix == 2)
        
        Vfr <- unlist(V)[seq(2, 4 * length(V), by = 4)]
        
        Ur <- U[1, ] * (random == 1) + U[2, ] * (random == 2)
        Uf <- U[1, ] * (fix == 1) + U[2, ] * (fix == 2)
        
        L <- sum(log(1 + Omega * Vrr)) + sum(Vrr/(1 + Omega * Vrr) * (muRandom - (Ur - psi * Vfr)/Vrr)^2 - (Ur - psi * Vfr)^2/Vrr - 2 * 
            psi * Uf + Vff * psi^2)
        
    }
    
    if (length(random) == 2) {
        print("wrong argument random")
    }
    return(L)
}




