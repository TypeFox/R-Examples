#'@title PANIC (2010) 
#'
#'@description This function performs the tests of PANIC (2010).
#  Each test has a critical value of -1.64 at a 5% significance level.
#'
#'@usage panic10(x, nfac, k1, jj, demean)
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
#'@return rho1 Estimation of the Pooled Autoregressive Coefficient.
#'
#'@return Model This function shows MP models A, B, and C. A assumes no deterministic component.
#' B assumes a constant and allows for a fixed effect. C allows a fixed effect
#' and trend.
#'
#'@return MP.tests If demean is false, a matrix containing results for models A and B.
#' If demean is true, returns results for Pa, Pb, and Model C
#'
#'@return test.C.P A matrix containing the Pooled tests as well as Model C.
#'
#'@return extra.test If demean is false, a matrix containing results PMSB test with fixed effects,
#'  rho, and the original ADF test on the idiosyncratic component from PANIC (2004).
#' If demean is true, returns results for PMSB test with fixed effects and trend, rho,
#'  and original LM test on idiosyncratic component from PANIC (2004)
#'
#'@references Bai, Jushan, and Serena Ng.
#''Panel Unit Root Tests With Cross-Section Dependence: A Further Investigation.'
#' Econometric Theory 26.04 (2010): 1088-1114. Print.

panic10 <- function(x, nfac, k1, jj, demean) {
    
    if (demean == FALSE) {
        x <- as.matrix(x)
        
        
        
        dx <- trimr(mydiff(x, 1), 1, 0)
        
        Tn <- dim(dx)[1]
        
        N <- dim(dx)[2]
        
        scale <- sqrt(N) * Tn
        
        factors <- getnfac(dx, nfac, jj)
        
        ic <- factors$ic
        
        if (ic == 0) ic <- 1 
        
        lamhat <- factors$lambda
        
        dfhat <- factors$Fhat
        
        fhat <- apply(dfhat, 2, cumsum)
        
        if (sum(lamhat) == 0) {
            dehat <- dx
        } else {
            dehat <- dx - dfhat %*% t(lamhat)
        }
        ehat0 <- apply(dehat, 2, cumsum)
        
        lagehat0 <- trimr(lagn(ehat0, 1), 1, 0)
        
        ehat0 <- trimr(ehat0, 1, 0)
        
        # Do old panic
        
        adf30 <- adf(ehat0, k1, -1)
        
        Pool <- pool(adfnc, t(adf30))
        
        adf30a <- Pool$adf31a
        
        adf30b <- Pool$adf31b
        
        # Model A compute rho0 (no demeaning)
        
        top0 <- sum(sum(lagehat0 * ehat0))
        
        bottom0 <- sum(sum(lagehat0 * lagehat0))
        
        rho0 <- top0/bottom0
        
        res0 <- ehat0 - lagehat0 * rho0
        
        Nuisance <- nuisance(res0, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        
        # tests using rho- (do not project on deterministic trends)
        
        ADJ <- N * Tn * HALF
        
        A1 <- 2
        
        B1 <- 1
        
        U1 <- (1/2)
        
        V1 <- (1/3)
        
        
        rho1 <- (top0 - ADJ)/bottom0
        # P = 0, -1 MP Tests Model A
        t_a <- scale * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
        
        t_b <- scale * (rho1 - 1) * sqrt(bottom0/(scale^2)) * sqrt(B1 * OMEGA2/PHI4)
        # P = 0 , -1 PMSB test
        t_c <- sqrt(N) * (sum(diag(ehat0 %*% t(ehat0)))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Model B tests that project on constant
        
        
        one <- matrix(1, I(Tn - 1), 1)
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        ehat <- Q_T %*% ehat0
        
        lagehat <- Q_T %*% lagehat0
        
        top <- sum(sum(lagehat * ehat))
        
        bottom <- sum(sum(lagehat * lagehat))
        
        rho1 <- top/bottom
        
        res1 <- ehat - lagehat * rho1
        
        Nuisance <- nuisance(res1, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        A1 <- 3
        
        B1 <- 2
        
        ADJ <- -N * Tn * SIG2/2
        
        rho1 <- (top - ADJ)/bottom
        
        # Model B for P = 0
        t_a1 <- scale * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
        
        t_a2 <- scale * (rho1 - 1) * sqrt(bottom/(scale^2)) * sqrt(B1 * OMEGA2/PHI4)
        
        
        
        
        
        test.A.B <- data.frame(Model.A=c( t_a, t_b), Model.B=c( t_a1, t_a2))
        
        rownames(test.A.B) <- c( "ta", "tb")
        
        extra.test <- data.frame(PMSB=t_c,rho1= rho1,Pool.ADF= adf30b)
        
        output <- list(MP.tests = test.A.B, PAC.tests = extra.test)
        
        return(output)
    } else {
        
        
        x <- as.matrix(x)
        
        dX <- trimr(mydiff(x, 1), 1, 0)
        
        intdX <- as.matrix(t(apply(dX, 2, mean)))
        
        repmat <- intdX[rep(seq_len(nrow(intdX)), each = I(nrow(x) - 1)), ]
        
        dx <- dX - repmat
        
        Tn <- dim(dx)[1]
        
        N <- dim(dx)[2]
        
        scale <- sqrt(N) * Tn
        
        factors <- getnfac(dx, nfac, jj)
        
        ic <- factors$ic
        
        if (ic == 0) ic <- 1
        
        lamhat <- factors$lambda
        
        dfhat <- as.matrix(factors$Fhat)
        
        fhat <- apply(dfhat, 2, cumsum)
        
        if (sum(sum(lamhat)) == 0) {
            
            dehat <- dx
        } else {
            
            dehat <- dx - dfhat %*% t(lamhat)
        }
        
        ehat0 <- apply(dehat, 2, cumsum)
        
        lagehat0 <- trimr(lagn(ehat0, 1), 1, 0)
        
        ehat0 <- trimr(ehat0, 1, 0)
        
        # Do old panic
        
        adf31 <- adf(ehat0, k1, -1)
        
        Pool <- pool(lm1, t(adf31))
        
        adf31a <- Pool$adf31a
        
        adf31b <- Pool$adf31b
        
        bottom0 <- sum(sum(lagehat0 * lagehat0))
        
        top0 <- sum(sum(lagehat0 * ehat0))
        
        rho0 <- top0/bottom0
        
        res0 <- ehat0 - lagehat0 * rho0
        
        Nuisance <- nuisance(res0, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        
        # No longer do detrending
        
        
        
        ADJ <- SIG2/OMEGA2
        
        A1 <- 36/5
        
        B1 <- 5/6
        
        U1 <- 1/6
        
        V1 <- 1/45
        
        
        # P = 1 for Pa and Pb
        t_a <- scale * (rho0 - 1 + ADJ * 3/Tn)/sqrt(A1 * PHI4 * SIG2^2/(OMEGA2^4))
        
        t_b <- scale * (rho0 - 1 + ADJ * 3/Tn) * sqrt(bottom0/(scale^2)) * sqrt(B1 * (OMEGA2^3)/(PHI4 * (SIG2^2)))
        # P = 1 PMSB
        t_c <- sqrt(N) * (sum(diag(ehat0 %*% t(ehat0)))/(N * Tn^2) - U1 * OMEGA2)/sqrt(V1 * PHI4)
        
        # Tests that project on intercept and trends
        
        one <- cbind(matrix(1, I(Tn - 1), 1), as.matrix(seq(1, I(Tn - 1))))
        
        Q_T <- diag(I(Tn - 1)) - one %*% solve(crossprod(one)) %*% t(one)
        
        ehat <- Q_T %*% ehat0
        
        lagehat <- Q_T %*% lagehat0
        
        top <- sum(sum(lagehat * ehat))
        
        bottom <- sum(sum(lagehat * lagehat))
        
        rho1 <- (top)/bottom
        
        res1 <- ehat - lagehat * rho1
        
        Nuisance <- nuisance(res1, 0)
        
        sig2 <- Nuisance$sig2
        
        omega2 <- Nuisance$omega2
        
        half <- Nuisance$half
        
        OMEGA2 <- mean(omega2)
        
        PHI4 <- mean(omega2 * omega2)
        
        SIG2 <- mean(sig2)
        
        HALF <- mean(half)
        
        A1 <- 15/4
        
        B1 <- 4
        
        ADJ <- -N * Tn * SIG2/2
        
        rho1 <- (top - ADJ)/bottom
        # Model C
        t_a1 <- scale * (rho1 - 1)/sqrt(A1 * PHI4/(OMEGA2 * OMEGA2))
        
        t_a2 <- scale * (rho1 - 1) * sqrt(bottom/(scale^2)) * sqrt(B1 * OMEGA2/PHI4)
        
        
        
        test.C.P <- data.frame(Pooled=c("Pa", "Pb"),Pool.value=c(t_a, t_b),C=c("ta", "tb"),C.Value=c(  t_a1,   t_a2))
        
        
        extra.test <- data.frame(PMSB = t_c, rho1 = rho1, Pool.LM = adf31b)
        
        
        output <- list(MP.tests = test.C.P, PAC.tests = extra.test)
        
        return(output)
        
    }
} 
