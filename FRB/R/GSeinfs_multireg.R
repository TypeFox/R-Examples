`GSeinfs_multireg` <-
function(X,Y,ests=GSest_multireg(X, Y))
{

# empirical influences for GS-estimates of multivariate regression

#-------------------------------------------------------------
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

X <- as.matrix(X)
p<-ncol(X)
int<-FALSE
interceptdetection <- apply(X==1, 2, all)
if (any(interceptdetection)) int<-TRUE
zonderint <- (1:p)[interceptdetection==FALSE]
Xzonderint <- X[,zonderint,drop=FALSE]
X <- as.matrix(Xzonderint)

n <- nrow(X)
p <- ncol(X)
m <- ncol(Y)

c <- ests$c
b <- ests$b
beta <- ests$coefficients

if(nrow(beta)>p) beta <- beta[2:(p+1),,drop=FALSE]
Sigma <- ests$Sigma

resmat <- Y-X%*%beta

places <- t(combn(1:n,2))
term1resvector <- as.matrix(resmat[places[,1],])
term2resvector <- as.matrix(resmat[places[,2],])
diffresmat <- term1resvector-term2resvector

divec <- sqrt(mahalanobis(diffresmat,rep(0,m),Sigma))
#psivec <- rhobiweightder1(divec,c)
psidervec <- rhobiweightder2(divec,c)
uvec <- rhobiweightder1(divec,c)/divec
vvec <- rhobiweightder1(divec,c) * divec

betacst <- (1-1/m)*mean(uvec)+1/m*mean(psidervec)
gamma1 <- (mean(psidervec*(divec^2))+(m+1)*mean(vvec))/(m+2)
gamma3 <- mean(vvec)
sigmaxinv <- solve(crossprod(X)/n)

infmatbeta <- matrix(0,n,p*m)
infmatsigma <- matrix(0,n,m*m)
sigmaypart <- matrix(0,n-1,m*m)

for (i in 1:n) {
    resimat <- matrix(rep(resmat[i,],n-1),ncol=m,byrow=TRUE)
    resmatzonder <- resmat[-i,]
    verschilmat <- resmatzonder-resimat
    divectimat <- sqrt(mahalanobis(verschilmat,rep(0,m),Sigma))
    uveci <- rhobiweightder1(divectimat,c)/divectimat
    party <- (1/betacst)*crossprod(uveci,verschilmat)/(n-1)
    infmatbeta[i,] <- vecop(sigmaxinv%*%(-as.matrix(X[i,])%*%party))
    for (j in 1:(n-1))
        {sigmaypart[j,]<- vecop((2/gamma1)*m*uveci[j]*tcrossprod(verschilmat[j,])-((2/gamma1)*rhobiweightder1(divectimat[j],c)*divectimat[j]-(4/gamma3)*(rhobiweight(divectimat[j],c)-b))*Sigma)}
       
    infmatsigma[i,] <- colMeans(sigmaypart)

}

return(list(infbeta=infmatbeta,infsigma=infmatsigma))
}

