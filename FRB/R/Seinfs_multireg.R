Seinfs_multireg <- function(X, Y, ests=Sest_multireg(X, Y)) {

# empirical influences for S-estimates of multivariate regression

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

c0 <- ests$c
b <- ests$b
covEst <- ests$Sigma
BetaEst <- ests$coefficients

sigmaXinv <- solve(crossprod(X)/n)
Sres <- Y - X %*% BetaEst
divec <- sqrt(mahalanobis(Sres, rep(0,m), covEst))

psidervecS <- rhobiweightder2(divec, c0)
uvecS <- rhobiweightder1(divec, c0) / divec
vvecS <- rhobiweightder1(divec,c0) * divec
rhovecS <- rhobiweight(divec,c0)

betaS <- (1-1/m) * mean(uvecS) + 1/m * mean(psidervecS)

gamma3 <- mean(vvecS)
gamma1S <- mean(psidervecS * (divec^2) + (m+1) * vvecS) / (m+2)

einfsbetaS <- matrix(0,n,p*m)
einfscovS <- matrix(0,n,m*m)
for (i in 1:n) {
    IFpart <- sigmaXinv %*% tcrossprod(X[i,], Sres[i,])
    IFbetaS <- 1 / betaS * uvecS[i] * IFpart
    einfsbetaS[i,] <- vecop( IFbetaS )

    IFpartcov <- tcrossprod(Sres[i,]) / divec[i]^2 - 1/m * covEst
    IFcovS <- 2/gamma3 * (rhovecS[i] - b) * covEst + 1/gamma1S * m * vvecS[i] * IFpartcov
    einfscovS[i,] <- vecop( IFcovS )
}
return(list(BetaS=einfsbetaS, covS=einfscovS))
 
}
