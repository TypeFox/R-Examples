###
### Assumption of nonignorability
###

NInocov <- function(Y, D, Z = NULL, CI = 0.95) {
  
  if (is.null(Z)) {
    ## Without noncompliance
    n <- length(D)
    R <- (!is.na(Y))*1
    pi00 <- mean((D == 0)*(R == 0))
    pi01 <- mean((D == 0)*(R == 1))
    pi10 <- mean((D == 1)*(R == 0))
    pi11 <- mean((D == 1)*(R == 1))
    
    p01 <- mean(Y[D == 0 & R == 1])
    p11 <- mean(Y[D == 1 & R == 1])
    
    p00 <- (p01*((1-p11)*pi00*pi11-(1-p01)*pi01*pi10))/(pi00*pi11*(p01-p11))
    p10 <- (p11*((1-p11)*pi00*pi11-(1-p01)*pi01*pi10))/(pi01*pi10*(p01-p11))
    
    p0 <- (p00*pi00+p01*pi01)/(pi00+pi01)
    p1 <- (p10*pi10+p11*pi11)/(pi10+pi11)
    Ybar <- matrix(c(p0, p1), nrow = 1)
    colnames(Ybar) <- c("Y0bar", "Y1bar")
    
    ATE <- p1-p0
    
    pi111 <- mean(Y == 1 & D == 1 & R == 1)
    pi101 <- mean(Y == 1 & D == 0 & R == 1)

    ## bounds
    lb <- (p11*pi11*(pi00+pi01)-(pi00+p01*pi01)*(pi10+pi11))/((pi10+pi11)*(pi00+pi01))
    ub <- ((pi10+p11*pi11)*(pi00+pi01)-p01*pi01*(pi10+pi11))/((pi10+pi11)*(pi00+pi01))
    
    ## asymptotic variance of ATE
    Sigma <- matrix(NA, ncol = 5, nrow = 5)
    Sigma[1,1] <- pi01*(1-pi01)
    Sigma[2,2] <- pi10*(1-pi10)
    Sigma[3,3] <- pi11*(1-pi11)
    Sigma[4,4] <- pi11*p11*(1-pi11*p11)
    Sigma[5,5] <- pi01*p01*(1-pi01*p01)
    Sigma[1,2] <- Sigma[2,1] <- -pi01*pi10
    Sigma[1,3] <- Sigma[3,1] <- -pi01*pi11
    Sigma[1,4] <- Sigma[4,1] <- -pi01*pi11*p11
    Sigma[1,5] <- Sigma[5,1] <- pi01*p01*(1-pi01)
    Sigma[2,3] <- Sigma[3,2] <- -pi10*pi11
    Sigma[2,4] <- Sigma[4,2] <- -pi10*pi11*p11
    Sigma[2,5] <- Sigma[5,2] <- -pi10*p01*pi01
    Sigma[3,4] <- Sigma[4,3] <- pi11*p11*(1-pi11)
    Sigma[3,5] <- Sigma[5,3] <- -pi11*pi01*p01
    Sigma[4,5] <- Sigma[5,4] <- -pi11*p11*pi01*p01
    
    Delta <- matrix(NA, ncol = 1, nrow = 5)
    A <- pi11*p11*(pi00+pi01)-pi01*p01*(pi10+pi11)
    B <- pi00*pi11-pi10*pi01
    C <- pi11*pi01*(p01-p11)*(pi00+pi01)*(pi10+pi11)
    Delta[1] <- -(pi10+pi11)*A+A*(B-A)*p11/(pi01*(p01-p11))
    Delta[2] <- -(p11*pi11+p01*pi01)*(B-2*A)-(pi11+pi01)*A -
      A*(B-A)*(pi00+pi01-pi10-pi11)/((pi00+pi01)*(pi10+pi11))
    Delta[3] <- -(p11*pi11+p01*pi01)*(B-2*A)+(pi00-pi11)*A -
      A*(B-A)*(p01/((p01-p11)*pi11)+(pi00+pi01-pi10-pi11)/((pi00+pi01)*(pi10+pi11)))
    Delta[4] <- (pi00+pi01)*(B-2*A)+A*(B-A)/(pi11*(p01-p11))
    Delta[5] <- -(pi10+pi11)*(B-2*A)-A*(B-A)/(pi01*(p01-p11))
    Delta <- Delta/C
    
    var <- t(Delta)%*%Sigma%*%Delta
    se <- sqrt(var/n)
    
    ps <- matrix(c(p00, p10, p01, p11), nrow = 1)
    colnames(ps) <- c("p00", "p10", "p01", "p11")
    
    pis <- matrix(c(pi00, pi10, pi01, pi11), nrow = 1)
    colnames(pis) <- c("pi00", "pi10", "pi01", "pi11")
    
    ATE <- matrix(c(ATE, se), nrow = 1)
    colnames(ATE) <- c("ATE", "s.e.")
    z <- qnorm(1-(1-CI)/2)
    ci.est <- matrix(c(ATE[1]-z*se, ATE[1]+z*se), nrow = 1)
    colnames(ci.est) <- c("lower CI", "upper CI")
    return(list(Ybar = Ybar, ATE = ATE, CI = ci.est,
                bounds = c(lb, ub), ps = ps, pis = pis, n = n)) 
  } else {
    ## WITH NONCOMPLIANCE
    n <- length(D)
    R <- (!is.na(Y))*1
    pi000 <- mean(D == 0 & R == 0 & Z == 0)
    pi100 <- mean(D == 1 & R == 0 & Z == 0)
    pi001 <- mean(D == 0 & R == 0 & Z == 1)
    pi101 <- mean(D == 1 & R == 0 & Z == 1)
    pi010 <- mean(D == 0 & R == 1 & Z == 0)
    pi110 <- mean(D == 1 & R == 1 & Z == 0)
    pi011 <- mean(D == 0 & R == 1 & Z == 1)
    pi111 <- mean(D == 1 & R == 1 & Z == 1)
    
    p010 <- mean(Y[D == 0 & R == 1 & Z == 0])
    p011 <- mean(Y[D == 0 & R == 1 & Z == 1])
    p110 <- mean(Y[D == 1 & R == 1 & Z == 0])
    p111 <- mean(Y[D == 1 & R == 1 & Z == 1])
    
    p000 <- (p010*((1-p011)*pi000*pi011-(1-p010)*pi010*pi001))/(pi000*pi011*(p010-p011))
    p100 <- (p110*((1-p111)*pi100*pi111-(1-p110)*pi110*pi101))/(pi100*pi111*(p110-p111))
    p001 <- (p011*((1-p011)*pi000*pi011-(1-p010)*pi010*pi001))/(pi010*pi001*(p010-p011))
    p101 <- (p111*((1-p111)*pi100*pi111-(1-p110)*pi110*pi101))/(pi110*pi101*(p110-p111))
    
    Y0bar <- (p000*pi000+p010*pi010+p100*pi100+p110*pi110)/(pi000+pi010+pi100+pi110)
    Y1bar <- (p001*pi001+p011*pi011+p101*pi101+p111*pi111)/(pi001+pi011+pi101+pi111)
    Ybar <- matrix(c(Y0bar, Y1bar), nrow = 1)
    colnames(Ybar) <- c("Y0bar", "Y1bar")
    
    ITT <- Y1bar-Y0bar
    CACE <- ITT/(mean(D[Z==1])-mean(D[Z==0]))
    
    ## alternative parameterization
    A0 <- p011*(pi011+((1-p011)*pi000*pi011/pi010-(1-p010)*pi001)/(p010-p011)) 
    A1 <- p111*(pi111+((1-p111)*pi100*pi111/pi110-(1-p110)*pi101)/(p110-p111))
    A <- A0 + A1
    C1 <- pi001+pi011+pi101+pi111
    C0 <- pi000+pi010+pi100+pi110
    Y1bar.alt <- A/C1
    
    B0 <- p010*(pi010+((1-p011)*pi000-(1-p010)*pi010*pi001/pi011)/(p010-p011)) 
    B1 <- p110*(pi110+((1-p111)*pi100-(1-p110)*pi110*pi101/pi111)/(p110-p111)) 
    B <- B0 + B1
    Y0bar.alt <- B/C0

    cace <- (Y1bar.alt-Y0bar.alt)*C0*C1
    cace.denom <- (pi101+pi111)*C0-(pi100+pi110)*C1
    cace <- cace/cace.denom

    ## variances
    Sigma <- matrix(NA, ncol = 11, nrow = 11)
    Sigma[1,1] <- pi001*(1-pi001)
    Sigma[2,2] <- pi010*(1-pi010)
    Sigma[3,3] <- pi100*(1-pi100)
    Sigma[4,4] <- pi011*(1-pi011)
    Sigma[5,5] <- pi101*(1-pi101)
    Sigma[6,6] <- pi110*(1-pi110)
    Sigma[7,7] <- pi111*(1-pi111)
    Sigma[8,8] <- p010*pi010*(1-p010*pi010)
    Sigma[9,9] <- p110*pi110*(1-p110*pi110)
    Sigma[10,10] <- p011*pi011*(1-p011*pi011)
    Sigma[11,11] <- p111*pi111*(1-p111*pi111)
    
    ## covariances
    Sigma[1,2] <- Sigma[2,1] <- -pi001*pi010
    Sigma[1,3] <- Sigma[3,1] <- -pi001*pi100
    Sigma[1,4] <- Sigma[4,1] <- -pi001*pi011
    Sigma[1,5] <- Sigma[5,1] <- -pi001*pi101
    Sigma[1,6] <- Sigma[6,1] <- -pi001*pi110
    Sigma[1,7] <- Sigma[7,1] <- -pi001*pi111
    Sigma[1,8] <- Sigma[8,1] <- -pi001*p010*pi010
    Sigma[1,9] <- Sigma[9,1] <- -pi001*p110*pi110
    Sigma[1,10] <- Sigma[10,1] <- -pi001*p011*pi011
    Sigma[1,11] <- Sigma[11,1] <- -pi001*p111*pi111
    
    Sigma[2,3] <- Sigma[3,2] <- -pi010*pi100
    Sigma[2,4] <- Sigma[4,2] <- -pi010*pi011
    Sigma[2,5] <- Sigma[5,2] <- -pi010*pi101
    Sigma[2,6] <- Sigma[6,2] <- -pi010*pi110
    Sigma[2,7] <- Sigma[7,2] <- -pi010*pi111
    Sigma[2,8] <- Sigma[8,2] <- p010*pi010*(1-pi010)
    Sigma[2,9] <- Sigma[9,2] <- -pi010*p110*pi110
    Sigma[2,10] <- Sigma[10,2] <- -pi010*p011*pi011
    Sigma[2,11] <- Sigma[11,2] <- -pi010*p111*pi111
    
    Sigma[3,4] <- Sigma[4,3] <- -pi100*pi011
    Sigma[3,5] <- Sigma[5,3] <- -pi100*pi101
    Sigma[3,6] <- Sigma[6,3] <- -pi100*pi110
    Sigma[3,7] <- Sigma[7,3] <- -pi100*pi111
    Sigma[3,8] <- Sigma[8,3] <- -pi100*p010*pi010
    Sigma[3,9] <- Sigma[9,3] <- -pi100*p110*pi110
    Sigma[3,10] <- Sigma[10,3] <- -pi100*p011*pi011
    Sigma[3,11] <- Sigma[11,3] <- -pi100*p111*pi111
    
    Sigma[4,5] <- Sigma[5,4] <- -pi011*pi101
    Sigma[4,6] <- Sigma[6,4] <- -pi011*pi110
    Sigma[4,7] <- Sigma[7,4] <- -pi011*pi111
    Sigma[4,8] <- Sigma[8,4] <- -pi011*p010*pi010
    Sigma[4,9] <- Sigma[9,4] <- -pi011*p110*pi110
    Sigma[4,10] <- Sigma[10,4] <- p011*pi011*(1-pi011)
    Sigma[4,11] <- Sigma[11,4] <- -pi011*p111*pi111
    
    Sigma[5,6] <- Sigma[6,5] <- -pi101*pi110
    Sigma[5,7] <- Sigma[7,5] <- -pi101*pi111
    Sigma[5,8] <- Sigma[8,5] <- -pi101*p010*pi010
    Sigma[5,9] <- Sigma[9,5] <- -pi101*p110*pi110
    Sigma[5,10] <- Sigma[10,5] <- -pi101*p011*pi011
    Sigma[5,11] <- Sigma[11,5] <- -pi101*p111*pi111
    
    Sigma[6,7] <- Sigma[7,6] <- -pi110*pi111
    Sigma[6,8] <- Sigma[8,6] <- -pi110*p010*pi010
    Sigma[6,9] <- Sigma[9,6] <- pi110*p110*(1-pi110)
    Sigma[6,10] <- Sigma[10,6] <- -pi110*p011*pi011
    Sigma[6,11] <- Sigma[11,6] <- -pi110*p111*pi111
    
    Sigma[7,8] <- Sigma[8,7] <- -pi111*p010*pi010
    Sigma[7,9] <- Sigma[9,7] <- -pi111*p110*pi110
    Sigma[7,10] <- Sigma[10,7] <- -pi111*p011*pi011
    Sigma[7,11] <- Sigma[11,7] <- p111*pi111*(1-pi111)
    
    Sigma[8,9] <- Sigma[9,8] <- -p010*pi010*p110*pi110
    Sigma[8,10] <- Sigma[10,8] <- -p010*pi010*p011*pi011
    Sigma[8,11] <- Sigma[11,8] <- -p010*pi010*p111*pi111
    
    Sigma[9,10] <- Sigma[10,9] <- -p110*pi110*p011*pi011
    Sigma[9,11] <- Sigma[11,9] <- -p110*pi110*p111*pi111

    Sigma[10,11] <- Sigma[11,10] <- -p011*pi011*p111*pi111
    
    ## partial derivatives
    dA <- dB <- matrix(NA, ncol = 1, nrow = 11)
    pi1010 <- p010*pi010
    pi1110 <- p110*pi110
    pi1011 <- p011*pi011
    pi1111 <- p111*pi111
    
    A0denom <- pi011*pi1010/pi1011-pi010
    A1denom <- pi111*pi1110/pi1111-pi110
    numer0 <- (pi011-pi1011)*pi000-pi001*(pi010-pi1010)
    numer1 <- (pi111-pi1111)*pi100-pi101*(pi110-pi1110)
    dA[1] <- -(pi010-pi1010)/A0denom
    dA[2] <- (-pi001*A0denom+numer0)/(A0denom^2)
    dA[3] <- (pi111-pi1111)/A1denom
    dA[4] <- (A0denom*pi000-pi1010*numer0/pi1011)/(A0denom^2)
    dA[5] <- -(pi110-pi1110)/A1denom
    dA[6] <- (-pi101*A1denom+numer1)/(A1denom^2)
    dA[7] <- (pi100*A1denom-pi1110*numer1/pi1111)/(A1denom^2)
    dA[8] <- (pi001*A0denom-pi011*numer0/pi1011)/(A0denom^2)
    dA[9] <- (pi101*A1denom-pi111*numer1/pi1111)/(A1denom^2)
    dA[10] <- 1 + (-pi000*A0denom+pi011*pi1010*numer0/(pi1011^2))/(A0denom^2)
    dA[11] <- 1 + (-pi100*A1denom+pi111*pi1110*numer1/(pi1111^2))/(A1denom^2) 
    
    B0denom <- pi011-pi1011*pi010/pi1010
    B1denom <- pi111-pi1111*pi110/pi1110
    dB[1] <- -(pi010-pi1010)/B0denom
    dB[2] <- (-pi001*B0denom+pi1011*numer0/pi1010)/(B0denom^2)
    dB[3] <- (pi111-pi1111)/B1denom
    dB[4] <- (pi000*B0denom-numer0)/(B0denom^2)
    dB[5] <- -(pi110-pi1110)/B1denom
    dB[6] <- (-pi101*B1denom+pi1111*numer1/pi1110)/(B1denom^2)
    dB[7] <- (pi100*B1denom-numer1)/(B1denom^2)
    dB[8] <- 1 + (pi001*B0denom-pi1011*pi010*numer0/(pi1010^2))/(B0denom^2)
    dB[9] <- 1 + (pi101*B1denom-pi1111*pi110*numer1/(pi1110^2))/(B1denom^2)
    dB[10] <- (-pi000*B0denom+pi010*numer0/pi1010)/(B0denom^2)
    dB[11] <- (-pi100*B1denom+pi110*numer1/pi1110)/(B1denom^2)
    
    dC0 <- dC1 <- matrix(0, ncol = 1, nrow = 11)
    dC0[2] <- dC0[3] <- dC0[6] <- dC1[1] <- dC1[4] <- dC1[5] <- dC1[7] <- 1
    
    ## put together
    partial <- (dA*C1 - dC1*A)/(C1^2)-(dB*C0-dC0*B)/(C0^2)
    ITT.se <- sqrt(t(partial)%*%Sigma%*%partial/n)
    z <- qnorm(1-(1-CI)/2)
    ITT.ci <- matrix(c(ITT-z*ITT.se, ITT+z*ITT.se), nrow = 1)
    colnames(ITT.ci) <- c("lower CI", "upper CI")
    ITT <- matrix(c(ITT, ITT.se), nrow = 1)
    colnames(ITT) <- c("ITT est.", "ITT s.e.")
    
    ## variance for CACE
    partial.cace <- (partial*C1*C0 + ITT[1]*(dC1*C0+C1*dC0))*cace.denom
    partial.cace[c(1, 2, 4, 8:11)] <- partial.cace[c(1, 2, 4, 8:11)] -
      (dC0[c(1, 2, 4, 8:11)]*(pi101+pi111) -
       (pi100+pi110)*dC1[c(1, 2, 4, 8:11)])*ITT[1]*C1*C0
    partial.cace[c(3, 6)] <- partial.cace[c(3,6)] -
      (dC0[c(3,6)]*(pi101+pi111) - C1 - (pi100+pi110)*dC1[c(3,6)])*ITT[1]*C1*C0
    partial.cace[c(5, 7)] <- partial.cace[c(5,7)] -
      (C0 + dC0[c(5, 7)]*(pi101+pi111) - (pi100+pi110)*dC1[c(5,7)])*ITT[1]*C1*C0
    partial.cace <- partial.cace/(cace.denom^2)
    CACE.se <- sqrt(t(partial.cace)%*%Sigma%*%partial.cace/n)
    CACE.ci <- matrix(c(CACE-z*CACE.se, CACE+z*CACE.se), nrow = 1)
    colnames(CACE.ci) <- c("lower CI", "upper CI")
    CACE <- matrix(c(CACE, CACE.se), nrow = 1)
    colnames(CACE) <- c("CACE est.", "CACE s.e.")
    
    ## output
    ps <- matrix(c(p000, p100, p001, p101, p010, p011, p110, p111), nrow = 1)
    colnames(ps) <- c("p000", "p100", "p001", "p101", "p010", "p011",
                      "p110", "p111")
    
    pis <- matrix(c(pi000, pi100, pi001, pi101, pi010, pi110, pi011, pi111), nrow = 1)
    colnames(pis) <- c("pi000", "pi100", "pi001", "pi101", "pi010",
                       "pi110", "pi011", "pi111")
    
    return(list(Ybar = Ybar, ITT = ITT, CACE = CACE, ITT.ci = ITT.ci,
                CACE.ci = CACE.ci, ps = ps, pis = pis, n = n))
  }
}
