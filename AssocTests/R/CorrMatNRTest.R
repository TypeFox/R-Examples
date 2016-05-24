CorrMatNRTest <- function(y, g)
{
    n0 <- sum(g==0)
    n1 <- sum(g==1)
    n2 <- sum(g==2)
    n  <- n0 + n1 + n2

    Y0 <- y[g==0]
    Y1 <- y[g==1]
    Y2 <- y[g==2]

 # genotype 1 relative to genotype 0
    d0  <- 1:n0
    d1  <- (n0+1):(n0+n1)
    RK1  <- rank(c(Y0,Y1))
    f01 <- ( sum(RK1[d1]) - n1*(n1+1)/2 ) / (n0*n1)

 #genotype 2 relative to genotype 0
    d2  <- (n0+1):(n0+n2)
    RK2  <- rank(c(Y0,Y2))
    f02 <- ( sum(RK2[d2])- n2*(n2+1)/2 ) / (n0*n2)

# genotype 2 relative to genotype 1
    d3 <- (n1+1):(n1+n2)
    RK12 <- rank(c(Y1,Y2))
    f12 <- (sum(RK12[d3])-n2*(n2+1)/2 ) / (n1*n2)

# standard deviation
    temp1 <- rep(NA, n0)
    temp2 <- rep(NA, n0)
    for (i in 1:n0)
    {
        temp1[i] <- ( sum(Y0[i]<Y1)/n1 - 1/2 )^2
        temp2[i] <- ( sum(Y0[i]<Y2)/n2 - 1/2 )^2
    }
    temp3 <- rep(NA, n1)
    temp4 <- rep(NA, n2)
    for (j in 1:n1)
    {
        temp3[j] <- ( sum(Y0<Y1[j])/n0 - 1/2 )^2
    }
    for (k in 1:n2)
    {
        temp4[k] <- ( sum(Y0<Y2[k])/n0 - 1/2 )^2
    }
    temp5 <- rep(NA, n1)
    temp6 <- rep(NA, n1)
    for (j in 1:n1)
    {
        temp5[j] <- sum(Y0<Y1[j])/n0 - 1/2
        temp6[j] <- sum(Y1[j]<Y2)/n2 - 1/2
    }
    temp7 <- temp6^2
    temp8 <- rep(NA, n2)
    for(k in 1:n2)
    {
        temp8[k] <- (sum(Y1<Y2[k])/n1 - 1/2)^2
    }

#### the necessary variance expressions
    f01.var <- 1/(n0*n1) * ( (n1-1)/n0 *sum(temp1) + (n0-1)/n1*sum(temp3) + 1/4 )
    f02.var <- 1/(n0*n2) * ( (n2-1)/n0 *sum(temp2) + (n0-1)/n2*sum(temp4) + 1/4 )
    f12.var <- 1/(n1*n2) * ( (n2-1)/n1 *sum(temp7) + (n1-1)/n2*sum(temp8) + 1/4 )
    cov01.12   <- (1/n1)^2*sum(temp5*temp6)
    a1 <-  sqrt((n0+n1)/f01.var)
    a2 <-  sqrt((n1+n2)/f12.var)
    w1 <-  a1/(a1+a2)
    w2 <-  a2/(a1+a2)
    sigma.A <- w1^2*f01.var+ w2^2*f12.var+  2*w1*w2*cov01.12
    temp9 <- rep(NA, n2)
    for(k in 1:n2)
    {
         temp9[k] <- (sum(Y0<Y2[k])/n0 -1/2 )*(sum(Y1<Y2[k])/n1 - 1/2)
    }
    cov02.12 <- (1/n2)^2*sum(temp9)
    sigma.R <- (n0/(n0+n1))^2* f02.var + (n1/(n0+n1))^2* f12.var +2*n0*n1/((n0+n1)^2)*cov02.12
    temp10 <- rep(NA, n0)
    for(i in 1:n0)
    {
         temp10[i] <- (sum(Y0[i]<Y1)/n1 -1/2 )*(sum(Y0[i]<Y2)/n2 -1/2)
    }
    cov01.02 <- 1/(n0^2)*sum(temp10)
    sigma.D <- (n1/(n1+n2))^2*f01.var + (n2/(n1+n2))^2*f02.var + 2*n1*n2/((n1+n2)^2)*cov01.02

############## estimate of the variance between the NRT test statistics Z_R, Z_A, Z_D
### 1\ the covariance between Z_R, Z_A
    RA.term1 <- n0/(n0+n1)*w1*cov01.02
    RA.term2 <- n0/(n0+n1)*w2*cov02.12
    RA.term3 <- n1/(n0+n1)*w1*cov01.12
    RA.term4 <- n1/(n0+n1)*w2*f12.var
    cov.RA <- RA.term1 + RA.term2 + RA.term3 + RA.term4

  ### 2\ the covariance between Z_A, Z_D
    AD.term1 <- n1/(n1+n2)*w1*f01.var
    AD.term2 <- n2/(n1+n2)*w1*cov01.02
    AD.term3 <- n1/(n1+n2)*w2*cov01.12
    AD.term4 <- n2/(n1+n2)*w2*cov02.12
    cov.AD <- AD.term1 + AD.term2 + AD.term3 + AD.term4

  ### 3\ the covariance between Z_R, Z_D
    RD.term1 <- n0*n1/((n0+n1)*(n1+n2))*cov01.02
    RD.term2 <- n0*n2/((n0+n1)*(n1+n2))*f02.var
    RD.term3 <- n1*n1/((n0+n1)*(n1+n2))*cov01.12
    RD.term4 <- n1*n2/((n0+n1)*(n1+n2))*cov02.12
    cov.RD <- RD.term1 +  RD.term2 + RD.term3 + RD.term4

 ### the covariance matrix
    cor.RA <- cov.RA/sqrt(sigma.R*sigma.A)
    cor.AD <- cov.AD/sqrt(sigma.A*sigma.D)
    cor.RD <- cov.RD/sqrt(sigma.R*sigma.D)

    rho <- matrix(0, nrow=3, ncol=3)
    rho[1,2:3] <- c(cor.RA, cor.RD)
    rho[2,3] <- cor.AD
    cov.mat <- rho + t(rho)
    cov.mat <- cov.mat + diag(3)

    cov.mat
}
