#'@title MCMC PANIC (2010) Power Test
#'
#'@description This function performs the tests of PANIC (2010) with the common and idiosyncratic components
#'  being resampled with the mcmcfactanal() function from MCMCpack.
#'  
#'@usage MCMCpanic10(x, nfac, k1, jj, demean = FALSE, burn = 1000, mcmc = 10000, thin = 10,
#' verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0,
#'   a0 = 0.001, b0 = 0.001, std.var = TRUE)
#'
#'
#'@param x A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test.
#'
#'@param jj an Integer 1 through 8. Choices 1 through 7 are respectively, IC(1),
#' IC(2), IC(3), AIC(1), BIC(1), AIC(3), and BIC(3), respectively. Choosing 8
#' makes the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#'@param demean logical argument. If TRUE, function performs tests on demeaned
#' data. If FALSE, uses non-demeanded data generating process.
#'
#'@param burn The number of burn in iterators for the sampler
#'
#'@param mcmc The number of iterations in the sampler
#'
#'@param thin The thinning interval used in the simulation. mcmc must be divisible by this value.
#'
#'@param verbose A positive integer which determines whether or not the progress of the
#' sampler is printed to the screen. If verbose is greater than 0 the iteration 
#' number and the factor loadings and uniqueness are printed to the screen 
#' every verboseth iteration.
#'
#'@param seed The seed for the random number generator.
#'
#'@param lambda.start Starting values for the factor loading matrix Lambda.
#'
#'@param psi.start Starting values for the uniqueness
#'
#'@param l0 The means of the independent Normal prior on the factor loadings
#'
#'@param L0 A scalar or a matrix with the same dimensions as lambda. The precision (inverse variances)
#' of the independent Normal prior on the factor loadings.
#' 
#'@param a0 scalar or a k-vector. Controls the shape of the inverse Gamma prior on the uniqueness.
#'
#'@param b0 Controls the scale of the inverse Gamma prior on the uniqueness.
#'
#'@param std.var if TRUE the variables are rescaled to have zero mean and unit variance.
#' Otherwise, the variables are rescaled to have zero mean, but retain their observed variances
#' 
#'@return adf.mcmc A list of the MCMC samples of the test statistics. If demeaned is set to TRUE, adf.mcmc
#' will have the tests Pa, Pb, Model C, PMSB, and rho1. If FALSE, adf.mcmc will have Model A, Model B, 
#' PMSB, and rho. All tests are degenerative and have a critical value of -1.64. 
#'
#'@references Bai, Jushan, and Serena Ng.
#''Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation.'
#' Econometric Theory 26.04 (2010): 1088-1114. Print.
#' 
#'@references Andrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'


MCMCpanic10 <- function(x, nfac, k1, jj, demean = FALSE, burn = 1000, mcmc = 10000, thin = 10, verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
    a0 = 0.001, b0 = 0.001, std.var = TRUE) {
    
    
    if (demean == FALSE) {
        x <- as.matrix(x)
        
        
        
        dx <- trimr(mydiff(x, 1), 1, 0)
        
        Tn <- dim(dx)[1]
        
        N <- dim(dx)[2]
        
        scale <- sqrt(N) * Tn
        
        factors <- getnfac(dx, nfac, jj)
        
        ic <- factors$ic
        
        if (ic == 0) ic <- 1
        
        fac.test <- MCMCfactanal(~., factors = ic, data = as.data.frame(dx), burnin = burn, mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, lambda.start = lambda.start, 
            psi.start = psi.start, l0 = l0, L0 = L0, a0 = a0, b0 = b0, store.scores = TRUE, std.var = std.var)
       
        
        A1 <- 2
        
        B1 <- 1
        
        U1 <- (1/2)
        
        V1 <- (1/3)
 
        one <- matrix(1, I(Tn - 1), 1)
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        
        A2 <- 3
        
        B2 <- 2

        lamhat <- NULL
        dfhat  <- NULL
        top    <- NULL
        bottom <- NULL
        Nuisance <- NULL
        sig2   <- NULL
        omega2 <- NULL
        PHI4   <- NULL
        SIG2   <- NULL
        HALF   <- NULL
        ADJ    <- NULL
        rho1   <- NULL
        t_a    <- NULL
        t_b    <- NULL
        t_c    <- NULL
        t_a1   <- NULL
        t_a2   <- NULL

        for (i in 1:I(mcmc/thin)) {

            lamhat <- matrix(fac.test[i, 1:I(N * ic)], N, ic)
            
            
            
            dfhat <- matrix(fac.test[i, I(N * ic + N + 1):I((Tn) * ic + N * ic + N)], I(Tn), ic, byrow = TRUE)
            

            top <- sum(sum(trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0) * trimr(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1, 0)))
            
            bottom <- sum(sum(trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0) * trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1), 1, 0)))
            
            
            
            Nuisance <- nuisance(trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0) - trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1), 1, 0) * (top/bottom), 0)
            
            
            
            sig2 <- Nuisance$sig2
            
            omega2 <- Nuisance$omega2
            
            half <- Nuisance$half
            
            OMEGA2 <- mean(omega2)
            
            PHI4 <- mean(omega2 * omega2)
            
            SIG2 <- mean(sig2)
            
            HALF <- mean(half)
        
        
        # tests using rho- (do not project on deterministic trends)

            ADJ  <- N * Tn * HALF
            
            rho1[[i]] <- (top - ADJ)/bottom
            # P = 0, -1 MP Tests Model A
            t_a[[i]] <- scale * (rho1[[i]] - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
            
            t_b[[i]] <- scale * (rho1[[i]] - 1) * sqrt(bottom/(scale^2)) * sqrt(B1 * OMEGA2/PHI4)
            # P = 0 , -1 PMSB test
            t_c[[i]] <- sqrt(N) * (sum(diag(trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0) %*% t(trimr(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1, 0))))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Model B tests that project on constant
        

            top <- sum(sum((Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0)) * (Q_T %*% trimr(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1, 0))))
            
            bottom <- sum(sum((Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0)) * (Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1), 1, 0))))
            
            Nuisance <- nuisance((Q_T %*% trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0)) - Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0) * (top/bottom), 0)
            
            sig2 <- Nuisance$sig2
            
            omega2 <- Nuisance$omega2
            
            half <- Nuisance$half
            
            OMEGA2 <- mean(omega2)
            
            PHI4 <- mean(omega2 * omega2)
            
            SIG2 <- mean(sig2)
            
            HALF <- mean(half)
            
            ADJ <- -N * Tn * SIG2/2
            
            rho1[[i]] <- (top - ADJ)/bottom
            
            # Model B for P = 0
            t_a1[[i]] <- scale * (rho1[[i]] - 1)/sqrt(A2 * PHI4/(OMEGA2 * OMEGA2))
            
            t_a2[[i]] <- scale * (rho1[[i]] - 1) * sqrt(bottom/(scale^2)) * sqrt(B2 * OMEGA2/PHI4)
        

        }
        
        
        
        output <- cbind(t_a, t_b, t_a1, t_a2, t_c, rho1)
        
        colnames(output) <- c("model A ta", "model A tb", "model B ta", "model B tb", "PMSB", "rho")
        
        return(output)
        ######################################################### 
    } else {
        
      
        x <- as.matrix(x)
        
        dx <- trimr(mydiff(x, 1), 1, 0)
        
        intdX <- as.matrix(t(apply(dx, 2, mean)))
        
        repmat <- intdX[rep(seq_len(nrow(intdX)), each = I(nrow(x) - 1)), ]
        
        dx <- dx - repmat
        
        Tn <- dim(dx)[1]
        
        N <- dim(dx)[2]
        
        scale <- sqrt(N) * Tn
        
        factors <- getnfac(dx, nfac, jj)
        
        ic <- factors$ic
        
        if (ic == 0) ic <- 1
        fac.test <- MCMCfactanal(~., factors = ic, data = as.data.frame(dx), burnin = burn, mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, lambda.start = lambda.start, 
            psi.start = psi.start, l0 = l0, L0 = L0, a0 = a0, b0 = b0, store.scores = TRUE, std.var = std.var)
  
        
        A1 <- 36/5
        
        B1 <- 5/6
        
        U1 <- 1/6
        
        V1 <- 1/45
        
        A2 <- 15/4
        
        B2 <- 4
        one <- cbind(matrix(1, I(Tn - 1), 1), as.matrix(seq(1, I(Tn - 1))))
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        lamhat <- NULL
        dfhat  <- NULL
        top    <- NULL
        bottom <- NULL
        Nuisance <- NULL
        sig2   <- NULL
        omega2 <- NULL
        half   <- NULL
        OMEGA2 <- NULL
        PHI4   <- NULL
        SIG2   <- NULL
        HALF   <- NULL
        ADJ    <- NULL
        t_a    <- NULL
        t_b    <- NULL
        t_c    <- NULL
        t_a1   <- NULL
        t_a2   <- NULL
        adf31  <- NULL
        adf31b <- NULL
        rho0   <- NULL
        rho1   <- NULL
        for (i in 1:I(mcmc/thin)) {
            lamhat <- matrix(fac.test[i, 1:I(N * ic)], N, ic)
            
            
            
            dfhat <- matrix(fac.test[i, I(N * ic + N + 1):I((Tn) * ic + N * ic + N)], I(Tn), ic, byrow = TRUE)

          
           top <- sum(sum(trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0) * trimr(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1, 0)))
           
            bottom <- sum(sum(trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0) * trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1), 1, 0)))
            
            Nuisance <- nuisance(trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0) - trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 
                                                                                                                           1), 1, 0) * (top/bottom), 0)
            sig2 <- Nuisance$sig2
            
            omega2 <- Nuisance$omega2
            
            half <- Nuisance$half
            
            OMEGA2 <- mean(omega2)
            
            PHI4 <- mean(omega2 * omega2)
            
            SIG2 <- mean(sig2)
            
            HALF <- mean(half)
            
        
        # No longer do detrending

            ADJ <- SIG2/OMEGA2
            
            rho0 <- top/bottom
            # P = 1 for Pa and Pb
            t_a[[i]] <- scale * (rho0 - 1 + ADJ * 3/Tn)/sqrt(A1 * PHI4 * SIG2^2/(OMEGA2^4))
            
            t_b[[i]] <- scale * (rho0 - 1 + ADJ * 3/Tn) * sqrt(bottom/(scale^2)) * sqrt(B1 * (OMEGA2^3)/(PHI4 * (SIG2^2)))
            # P = 1 PMSB
            t_c[[i]] <- sqrt(N) * (sum(diag(trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0) %*% t(trimr(apply(dx - tcrossprod(dfhat, lamhat), 
                2, cumsum), 1, 0))))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Tests that project on intercept and trends

            top <- sum(sum((Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1),
                                               1, 0)) * (Q_T %*% trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1, 0))))
            
            bottom <- sum(sum((Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1),
                                                  1, 0)) * (Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0))))
            
            Nuisance <- nuisance(( Q_T %*% trimr(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum),
                                                      1, 0)) - ( Q_T %*% trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum),
                                                                                    1), 1, 0)) * (top/bottom), 0)
            sig2 <- Nuisance$sig2
            
            omega2 <- Nuisance$omega2
            
            half <- Nuisance$half
            
            OMEGA2 <- mean(omega2)
            
            PHI4 <- mean(omega2 * omega2)
            
            SIG2 <- mean(sig2)
            
            HALF <- mean(half)
  
            ADJ <- -N * Tn * SIG2/2
            
            rho1[[i]] <- (top - ADJ)/bottom
            # Model C
            t_a1[[i]] <- scale * (rho1[[i]] - 1)/sqrt(A2 * PHI4/(OMEGA2 * OMEGA2))
            
            t_a2[[i]] <- scale * (rho1[[i]] - 1) * sqrt(bottom/(scale^2)) * sqrt(B2 * OMEGA2/PHI4)
        
            adf31 <- adf(trimr(lagn(apply(dx - tcrossprod(dfhat, lamhat), 2, cumsum), 1), 1, 0), k1, -1)
        
            adf31b[[i]] <- pool(lm1, t(adf31))$adf31b
            
        }
        
        
        adf.mcmc <- cbind(t_a, t_b, t_a1, t_a2, t_c, rho1,adf31b)
        
        colnames(adf.mcmc) <- c("Pa", "Pb", "Model C ta", "Model C tb", "PMSB", "rho1","lm")
        
        output <- list(adf.mcmc = adf.mcmc)
        
        
        return(output)
        
    }
} 
