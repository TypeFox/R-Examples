#'@title PANIC (2004) Power Test
#'
#'@description This function performs an MCMC over PANIC (2010) 
#'
#'@usage MCMCpanic04(x, nfac, k1, jj,burn = 1000, mcmc = 10000, thin = 10, verbose = 0,
#'seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, 
#'a0 = 0.001, b0 = 0.001, std.var = TRUE)
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
#'@param burn The number of burn in iterators for the sampler
#'
#'@param mcmc The number of iterations in the sampler
#'
#'@param thin The thinning interval used in the simulation. mcmc must be divisible by this value.
#'
#'@param verbose A positive integer which determines whether or not the progress of the
#' sampler is printed to the screen. If verbose is greater than 0 the iteration 
#' number and the factor loadings and uniquenesses are printed to the screen 
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
#'@return adf.mcmc A list of the MCMC samples of the test statistics. Returns the test statistics
#'  Pooled Cointegration a,  Pooled Cointegration b, Pooled Idiosyncratic a, 
#'  Pooled Idiosyncratic b, Pooled Demeaned test, and tests on Common components. The critical values for the Pooled Cointegration test
#'  can be found in Bai and Ng (2004). The pooled idiosyncratic test has a critical
#'  value of 1.64. The Pooled Demeaned test has a critical value of 2.86. The common components have a critical
#'  value of -2.86. All critical values are at a 5% significance level
#'
#'@return factor_mcmc The MCMC object from MCMCfactanal()
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#' 
#' @references ndrew D. Martin, Kevin M. Quinn, Jong Hee Park (2011). MCMCpack: Markov Chain Monte Carlo
#' in R. Journal of Statistical Software. 42(9): 1-21. URL http://www.jstatsoft.org/v42/i09/.
#'
MCMCpanic04 <- function(x, nfac, k1, jj, burn = 1000, mcmc = 10000, thin = 10, verbose = 0, seed = NA, lambda.start = NA, psi.start = NA, l0 = 0, L0 = 0, a0 = 0.001, 
    b0 = 0.001, std.var = TRUE) {
   
  x <- as.matrix(x)
    
    Tn <- dim(x)[1]
    
    N <- dim(x)[2]
    
    intx <- as.matrix(t(apply(x, 2, mean)))
    
    repmat <- intx[rep(seq_len(nrow(intx)), each = I(nrow(x))), ]
    
    x1 <- x - repmat
    
    x2 <- x1[2:Tn, ]
    
    dx <- trimr(mydiff(x, 1), 1, 0)
    
    scale <- sqrt(N) * Tn
    
    factors <- getnfac(dx, nfac, jj)

    ic <- factors$ic
    
    if (ic == 0) {
        ic <- 1
    }
    
    PC <- pc(dx, ic)
    
    if (ic == 1) {
      p = 0
    }
    if (ic == 0) {
      p = -1
    }
    if (ic > 1) {
      p = 1
    }
    
    lamhat <- PC$lambda
    
    dfhat <- PC$fhat
    
    fac.test <- MCMCfactanal(~., factors = ic, data = as.data.frame(dx), burnin = burn, mcmc = mcmc, thin = thin, verbose = verbose, seed = seed, lambda.start = lambda.start, 
        psi.start = psi.start, l0 = l0, L0 = L0, a0 = a0, b0 = b0, store.scores = TRUE, std.var = std.var)
    
    

    lamhat <- NULL
    dfhat  <- NULL
    ehat0  <- NULL
    fhat0  <- NULL
    reg    <- NULL
    ehat1  <- NULL
    beta1  <- NULL
    adf20  <- NULL
    adf30  <- NULL
    adf50  <- NULL
    padf30 <- NULL
    adf30a <- NULL
    adf30b <- NULL
    padf50 <- NULL
    adf50a <- NULL
    adf50b <- NULL
    adf20ab <- NULL
  
    for (j in 1:I(mcmc/thin)) {
        lamhat[[j]] <- matrix(fac.test[j, 1:I(N * ic)], N, ic)
     
        dfhat[[j]] <- matrix(fac.test[j, I(N * ic + N + 1):I((Tn - 1) * ic + N * ic + N)], I(Tn - 1), ic, byrow = TRUE)

        ehat0[[j]] <- apply(dx - tcrossprod(dfhat[[j]], lamhat[[j]]), 2, cumsum)
        
        fhat0[[j]] <- apply(dfhat[[j]], 2, cumsum)
        
        reg[[j]] <- cbind(matrix(1, I(Tn - 1), 1), fhat0[[j]])
        
        ehat1[[j]] <- matrix(0, I(Tn - 1), N)
        
        beta1[[j]] <- matrix(0, I(ic + 1), N)
        
        for (i in 1:N) {
            
            beta1[[j]][, i] <- qr.solve(reg[[j]], x2[, i])
            
            ehat1[[j]][, i] <- x2[, i] - reg[[j]] %*% beta1[[j]][, i]
        }
        
 
        adf20[[j]] <- matrix(0, 1, ic)
        
        for (i in 1:ic) {
            
            adf20[[j]][i] <- adf04(fhat0[[j]][, i], k1, p)
        }
        
        
        adf30[[j]] <- adf04(ehat0[[j]], k1, -1)  # test ehat0
        
        adf50[[j]] <- adf04(ehat1[[j]], k1, -1)  # test ehat1
        
        # now do the pooled test
        
        
        padf30[[j]] <- pool(adfnc, adf30[[j]])
        
        adf30a[[j]] <- padf30[[j]]$adf31a
        
        adf30b[[j]] <- padf30[[j]]$adf31b
        
        padf50[[j]] <- poolcoint(coint0, adf50[[j]], ic)
        
        adf50a[[j]] <- padf50[[j]]$pvala
        
        adf50b[[j]] <- padf50[[j]]$pvalb
        
    }
    
    adf20ab <- as.data.frame(matrix(unlist(adf20), mcmc, ic, byrow = TRUE))
    
    for (i in 1:ic) {
        colnames(adf20ab)[i] <- paste0("Common", i)
    }
    
    adf.tests <- cbind(adf50a, adf50b, adf30a, adf30b, adf20ab)
    
    
    colnames(c("Pooled Cointegration a", "Pooled Cointegration b", "Pooled Idiosyncratic a", "Pooled Idiosyncratic b", "Pooled Demeaned"))
    results <- list(adf.mcmc = adf.tests, factor_MCMC = fac.test)
} 
