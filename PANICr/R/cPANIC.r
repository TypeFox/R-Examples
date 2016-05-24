#'@title PANIC (2004) Non-Stationarity Tests on Common and Idiosyncratic Components
#'
#'@description This function performs PANIC (2010) Model C, PAC, and PMSB tests.
#' PAC estimates the pooled autoregressive coefficient, PMSB uses a sample
#' moment, and Model C performs the MP test while projecting on intercept and trend.
#'  The sample moments test is based off of the modified Sargan-Bhargava test (PMSB).
#'
#'@usage panic04(x, nfac, k1, jj)
#'
#'
#'@param x A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@param k1 The maximum lag allowed in the ADF test. The suggested setting for the maximum number of lags is
#' 4*(T/100)^(1/4) where T is the number of observations.
#'
#'@param jj an Integer 1 through 8. Choices 1 through 7 are respectively, IC(1),
#' IC(2), IC(3), AIC(1), BIC(1), AIC(3), and BIC(3), respectively. Choosing 8
#' makes the number of factors equal to the number of columns whose sum of
#' eigenvalues is less than  or equal to .5.
#'
#'@return adff A data frame containing pooled demeaned critical values, demeaned error term critical values ,demeaned
#'  and detrended critical values ,R squared for principle component
#'  , and the significance of the error components.
#'
#'@return adf.ind A matrix containing the critical values for the pooled Demeaned ADF test on the data, the
#' pooled ADF test on the common components, the pooled demeaned ADF test on the
#' Idiosyncratic component, and the pooled first differenced and demeaned ADF test on the
#' Idiosyncratic component.
#' 
#'@references Bai, Jushan, and Serena Ng. 
#''A PANIC Attack on Unit Roots and Cointegration.'
#' Econometrica 72.4 (2004): 1127-1177. Print.
#'
panic04 <- function(x, nfac, k1, jj) {
    
    
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
    
    lamhat <- PC$lambda
    
    dfhat <- PC$fhat
    
    
    
    
    if (sum(sum(lamhat)) == 0) {
        
        dehat <- dx
    } else {
        
        dehat <- dx - tcrossprod(dfhat, lamhat)
    }
    
    
    
    fhat0 <- apply(dfhat, 2, cumsum)
    
    ehat0 <- apply(dehat, 2, cumsum)
    
    ehat1 <- matrix(0, I(Tn - 1), N)
    
    reg <- cbind(matrix(1, I(Tn - 1), 1), fhat0)
    
    beta1 <- matrix(0, I(ic + 1), N)
    
    for (i in 1:N) {
        
        beta1[, i] <- qr.solve(reg, x2[, i])
        
        ehat1[, i] <- x2[, i] - reg %*% beta1[, i]
    }
    
    
    # some diagnostics to see the importance of the factors
    
    R21 <- matrix(0, N, 1)
    
    R22 <- matrix(0, N, 1)
    
    fit <- matrix(0, dim(fhat0)[1], N)
    
    fit1 <- matrix(0, I(Tn - 1), N)
    
    fit2 <- matrix(0, I(Tn - 1), N)
    
    lamhat <- as.matrix(lamhat)
    for (i in 1:N) {
        
        fit1[, i] <- fhat0 %*% lamhat[i, ]
        
        fit2[, i] <- dfhat %*% lamhat[i, ]
        
        R21[i, ] <- sd(dehat[, i])^2/sd(dx[, i])^2
        
        R22[i, ] <- sd(fit1[, i])/sd(ehat0[, i])
    }
    
    
    
    
    adf10 <- matrix(0, N, 1)
    
    adf20 <- matrix(0, ic, 1)
    
    adf30 <- matrix(0, N, 1)
    
    adf40 <- matrix(0, N, 1)
    
    adf50 <- matrix(0, N, 1)
    
    if (ic == 1) {
        p = 0
    }
    if (ic == 0) {
        p = -1
    }
    if (ic > 1) {
        p = 1
    }
    p <- 1
  
    adf10 <- adf04(x1, k1, p)
    
    for (i in 1:ic) {
        
        adf20[i, ] <- adf04(fhat0[, i], k1, p)  # test fhat0 for a unit root
    }
    
    adf30 <- adf04(ehat0, k1, -1)  # test ehat0
    
    adf50 <- adf04(ehat1, k1, -1)  # test ehat1
    
    # now do the pooled test
    
    padf10 <- pool(adfc2, adf10)
    
    adf10a <- padf10$adf31a
    
    adf10b <- padf10$adf31b
    
    
    padf30 <- pool(adfnc, adf30)
    
    adf30a <- padf30$adf31a
    
    adf30b <- padf30$adf31b
    
    
    padf50 <- poolcoint(coint0, adf50, ic)
    
    adf50a <- padf50$pvala
    
    adf50b <- padf50$pvalb
    
    
    
    adfr <- as.data.frame(cbind(as.matrix(seq(1:N)), t(adf10), t(adf30), t(adf50), R21, R22))
    colnames(adfr) <- c("Series", "adf", "ehat", "ehat1", "R2", "sifF/sige")
    
    Common <- matrix(adf20)
    colnames(Common) <- c("Common Test")
    
    adf.ind <- data.frame(Test=c("Pooled Demeaned","Pooled Idiosyncratic","Pooled Cointegration test"),Value=c(adf10b,  adf30b , adf50b))

    results <- list(adff = adfr, pooladf = adf.ind, Common = Common)
    
    return(results)
    
} 
