#'@title Determining The Number of Factors In Approximate Factor Model
#'
#'@description This function approximates the number of factors in an approximate factor model
#'  for large N by T matrices using the methods found in Bai and Ng (2002)
#'
#'@usage getnfac(x,kmax,jj)
#'
#'@param x A NxT matrix containing the data.
#'
#'@param kmax The maximum number of common factors used to compute the criterion
#' function for the estimate of ic, the number of common factors. This methedology
#' is weak to underestimation of the number of common factors. It is suggested
#' that overestimation is preferred.
#'
#' @param jj an Integer 1 through 8. Choices 1 through 7 are respectively, IC(1),
#' IC(2), IC(3), AIC(1), BIC(1), AIC(3), and BIC(3), respectively. Choosing 8
#' makes the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#' @details This function approximates the number of factors in an approximate
#' factor model. Amongst the penalty functions BIC(3) has been found to be
#' strict against cross-sectional dependence and is recommended for large
#' matrices. IC(1), IC(2), and IC(3) . AIC(1) will not work for all N and T.
#' BIC(1) will not work for small N relative to T. AIC(3) and BIC(3) take into
#' account the panel structure of the data. AIC(3) performs in consistently
#' across configurations of the data while BIC(3) may not perform well for some
#' configurations.
#'
#' @return ic The approximate number of factors based off of the chosen.
#' penalty function
#'
#' @return lambda Estimated factor loadings associated with common factors.
#'
#' @return Fhat Estimated common component
#'
#'
#' @references Jushan, and Serena Ng. 'Determining the Number of Factors in
#' Approximate Factor Models.' Econometrica 70.1 (2002): 191-221. Print.
#'
#'
getnfac <- function(x, kmax, jj) {
    
    x <- as.matrix(x)
    Tn <- dim(x)[1]
    
    N <- dim(x)[2]
    
    NT <- N * Tn
    
    NT1 <- N + Tn
    
    CT <- matrix(0, 1, kmax)
    
    ii <- seq(1:kmax)
    
    GCT <- min(N, Tn)
    
    if (jj == 1) {
        CT[1, ] <- log(NT/NT1) * ii * NT1/NT
    }
    
    if (jj == 2) {
        CT[1, ] <- (NT1/NT) * log(GCT) * ii
    }
    
    if (jj == 3) {
        CT[1, ] <- ii * log(GCT)/GCT
    }
    
    if (jj == 4) {
        CT[1, ] <- 2 * ii/Tn
    }
    
    if (jj == 5) {
        CT[1, ] <- log(T) * ii/Tn
    }
    
    if (jj == 6) {
        CT[1, ] <- 2 * ii * NT1/NT
    }
    
    if (jj == 7) {
        CT[1, ] <- log(NT) * ii * NT1/NT
    }
    
    IC1 <- matrix(0, dim(CT)[1], I(kmax + 1))
    
    Sigma <- matrix(0, 1, I(kmax + 1))
    
    XX <- x %*% t(x)
    
    eig <- svd(t(XX))
    
    Fhat0 <- eig$u
    
    eigval <- as.matrix(eig$d)
    
    Fhat1 <- eig$v
    
    sumeigval <- apply(eigval, 2, cumsum)/sum(eigval)
    
    if (jj < 8) {
        for (i in kmax:1) {
            
            Fhat <- Fhat0[, 1:i]
            
            lambda <- crossprod(Fhat, x)
            
            chat <- Fhat %*% lambda
            
            ehat = x - chat
            
            Sigma[i] <- mean(sum(ehat * ehat/Tn))
            
            IC1[, i] <- log(Sigma[i]) + CT[, i]
        }
        
        Sigma[kmax + 1] <- mean(sum(x * x/Tn))
        
        IC1[, kmax + 1] <- log(Sigma[kmax + 1])
        
        ic1 <- minindc(t(IC1))
        
        ic1 <- ifelse(ic1 <= kmax, ic1 * 1, ic1 * 0)
        
    }
    
    if (jj == 8) {
        
        for (j in 1:I(nrow(sumeigval))) {
            
            if (sumeigval[j] >= 0.5) {
                ic1 = j
                break
            }
        }
        
    }
    
    if (ic1 == 0) {
        
        Fhat = matrix(0, T, N)
        
        lambda = matrix(0, N, T)
        
        chat = matrix(0, T, N)
    } else {
        
        Fhat <- Fhat0[, 1:ic1]
        
        lambda <- crossprod(x, Fhat)
        
        chat <- Fhat %*% t(lambda)
    }
    
    
    output <- list(ic = ic1, lambda = lambda, Fhat = Fhat)
    return(output)
} 
