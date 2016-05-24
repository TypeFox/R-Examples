Seinfs_pca <- function(Y, ests=Sest_loccov(Y)) {

# empirical influences for S-estimates + PCA estimates

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
p <- ncol(Y)

c0 <- ests$c
b <- ests$b

locEst <- ests$Mu
covEst <- ests$Sigma
GammaEst <- det(covEst)^(-1/p) * covEst

# compute eigenvalues/vectors
eigresGamma <- eigen(GammaEst)
IXGamma <- order(eigresGamma$values)
eigs <- eigresGamma$values[IXGamma]
eigvecs <- eigresGamma$vectors[,IXGamma]

for (k in 1:p) {
    IX  <- order(abs(eigvecs[,k]))
    eigvecs[,k] <- sign(eigvecs[IX[p],k]) * eigvecs[,k]
}    

varperc <- rep(0,p-1)
for (k in 1:(p-1)) {
    varperc[k] <- sum(eigs[(p-k+1):p]) / sum(eigs)
}

divec <- sqrt(mahalanobis(Y, locEst, covEst))

psidervecS <- rhobiweightder2(divec, c0)
uvecS <- rhobiweightder1(divec, c0) / divec
vvecS <- rhobiweightder1(divec,c0) * divec
rhovecS <- rhobiweight(divec,c0)

betaS <- (1-1/p) * mean(uvecS) + 1/p * mean(psidervecS)
einfslocS <- 1 / betaS * matrix(rep(uvecS,p), ncol=p) * (Y - matrix(rep(locEst,n), n, byrow=TRUE))

gamma3 <- mean(vvecS)
gamma1S <- mean(psidervecS * (divec^2) + (p+1) * vvecS) / (p+2)

einfscovS <- matrix(0,n,p*p)
einfsshapeS <- matrix(0,n,p*p)
einfseigvec <- matrix(0,n,p*p)
einfseigscov <- matrix(0,n,p)
einfseigs <- matrix(0,n,p)
einfsvarperc <- matrix(0,n,p-1)
for (i in 1:n) {
  
    IFcovS <- 2/gamma3 * (rhovecS[i] - b) * covEst + 1/gamma1S * p * vvecS[i] * (t(Y[i,]-locEst) %*% (Y[i,]-locEst) / divec[i]^2 - 1/p * covEst)
    einfscovS[i,] <- vecop( IFcovS )
    
    IFshapeS <- 1/gamma1S * p * vvecS[i] * det(covEst)^(-1/p) * (t(Y[i,]-locEst) %*% (Y[i,]-locEst) / divec[i]^2 - 1/p * covEst)
    einfsshapeS[i,] <- vecop( IFshapeS )
    
    IFeigvecs <- matrix(0,p,p)
    for (k in 1:p) { 
        IFeigveck <- rep(0,p)
        for (j in 1:p) {
            if (j != k) {
                IFeigveck <- IFeigveck + 1 / (eigs[k]-eigs[j]) * (t(eigvecs[,j]) %*% IFshapeS %*% eigvecs[,k]) %*% eigvecs[,j]
            }
        }
        IFeigvecs[,k] <- IFeigveck
        einfseigs[i,k] <- t(eigvecs[,k]) %*% IFshapeS %*% eigvecs[,k]
        einfseigscov[i,k] <- t(eigvecs[,k]) %*% IFcovS %*% eigvecs[,k]
    }
    einfseigvec[i,] <- vecop(IFeigvecs)
    for (k in 1:(p-1)) {
        einfsvarperc[i,k] <- 1 / sum(eigs) * ((1-varperc[k]) * sum(einfseigs[i,(p-k+1):p]) - varperc[k] * sum(einfseigs[i,1:(p-k)]))
    }
}
return(list(loc=einfslocS, shape=einfsshapeS, cov=einfscovS, eigvec=einfseigvec, eigs=einfseigs, eigscov=einfseigscov, varperc=einfsvarperc))
 
}
