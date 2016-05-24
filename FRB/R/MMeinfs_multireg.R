MMeinfs_multireg <- function(X, Y, ests = MMest_multireg(X, Y)) {

# empirical influences for MM-estimates

# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder1 <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder2 <- function(x,c)
{
# Computes derivative of Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 6*x^2/(c^2) + 5*x^4/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------
# -                         main function                            -
# --------------------------------------------------------------------

n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)

c0 <- ests$c0
b <- ests$b
c1 <- ests$c1
MMBeta <- ests$coefficients
MMSigma <- ests$Sigma

sigmaXinv <- solve(crossprod(X)/n)
MMres <- Y - X %*% MMBeta
divec <- sqrt(mahalanobis(MMres, rep(0,m), MMSigma))

psidervec <- rhobiweightder2(divec, c1)
psidervecS <- rhobiweightder2(divec, c0)
uvec <- rhobiweightder1(divec, c1) / divec
uvecS <- rhobiweightder1(divec, c0) / divec
vvec <- rhobiweightder1(divec,c1) * divec
vvecS <- rhobiweightder1(divec,c0) * divec
rhovecS <- rhobiweight(divec,c0)

betaMM <- (1-1/m) * mean(uvec) + 1/m * mean(psidervec)
betaS <- (1-1/m) * mean(uvecS) + 1/m * mean(psidervecS)

gamma3 <- mean(vvecS)
gamma1MM <- mean(psidervec * (divec^2) + (m+1) * vvec) / (m+2)
gamma1S <- mean(psidervecS * (divec^2) + (m+1) * vvecS) / (m+2)

einfsbeta <- matrix(0,n,p*m)
einfsbetaS <- matrix(0,n,p*m)
einfscov <- matrix(0,n,m*m)
einfsshape <- matrix(0,n,m*m)
einfscovS <- matrix(0,n,m*m)
for (i in 1:n) {
    IFpart <- sigmaXinv %*% tcrossprod(X[i,], MMres[i,])
    IFbeta <- 1 / betaMM * uvec[i] * IFpart
    einfsbeta[i,] <- vecop( IFbeta )
    
    IFbetaS <- 1 / betaS * uvecS[i] * IFpart
    einfsbetaS[i,] <- vecop( IFbetaS )

    IFpartcov <- tcrossprod(MMres[i,]) / divec[i]^2 - 1/m * MMSigma
    IFcov <- 2/gamma3 * (rhovecS[i] - b) * MMSigma + 1/gamma1MM * m * vvec[i] * IFpartcov     
    einfscov[i,] <- vecop( IFcov )
    
    IFcovS <- 2/gamma3 * (rhovecS[i] - b) * MMSigma + 1/gamma1S * m * vvecS[i] * IFpartcov
    einfscovS[i,] <- vecop( IFcovS )

    IFshape <- 1/gamma1MM * m * vvec[i] * det(MMSigma)^(-1/m) * IFpartcov
    einfsshape[i,] <- vecop( IFshape )
}
return(list(Beta=einfsbeta, BetaS=einfsbetaS, cov=einfscov, covS=einfscovS, shape=einfsshape))
 
}
