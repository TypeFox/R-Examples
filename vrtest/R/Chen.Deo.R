Chen.Deo <-
function(x,kvec)
{

ks <- max(kvec)
sig2 <- sd(x)^2
z <- as.matrix(x - mean(x))
n <- nrow(z)

lamda <- as.matrix(2*pi/n * (1:as.integer((n-1)/2)))


#@ equation 9 of Chen and Deo @
wmat <- matrix(NA,nrow=nrow(lamda),ncol=length(kvec))
for(i in 1:nrow(lamda)) 
{
wmat[i,] = (kvec^(-1) * ( sin( kvec*lamda[i]/2 ) /sin( lamda[i]/2 ) )^2 )
}
w1 <-colSums(wmat)
w2 <-colSums(wmat^2)
w3 <-colSums(wmat^3)

beta  <- 1 - 2/3*( w1*w3) /(w2^2)



#@ equation 10 of Chen and Deo @

tem <- complex(imaginary=-lamda)
Ilamda <- matrix(NA,nrow=nrow(lamda),1)
for( j in 1:nrow(lamda) )
{
    sum1 <- 0
    for( i in 1:n )
    {    sum1 <- sum1 + z[i] * exp(tem[j]*i) }
    
Ilamda[j] <- (2*pi*n)^(-1) * Mod(sum1)^2
}

sum1 <- numeric(0)
for( j in 1:length(kvec) )
{
sum1 <- c(sum1,colSums(Ilamda * wmat[,j]) * (1 - kvec[j]/n)^(-1) * (4*pi) / (n*sig2) )
}

Vp <- sum1^(beta)


#@ tau values on page 215 @
tauvec <- numeric(0)
for( j in 1:ks)
{
tauvec <- c(tauvec, sum( (z[(j+1):n]^2)*(z[1:(n-j)]^2))*sig2^(-2)/(n-j-4))
}


cnk <- n*(n-kvec)^(-1)

#@ Matrices in equation 11 @
lmat <- matrix(0,nrow=ks+1,ncol=length(kvec))

for( i in 1:length(kvec))
{
    for( j in 1:(kvec[i]-1))
    {lmat[j,i] <- 2*cnk[i]*(1-j/kvec[i])}
}
lmat[ks+1,] <- -(kvec * cnk - n/(n-1))

bvec <- matrix(0,nrow=ks,ncol=1)
for(j in 1:ks)
{bvec[j] <- 2*(n-j)*n^(-3)*tauvec[j] + 2*j*n^(-3)}

avec <- ( ( n- 1:ks ) * tauvec )/n^2 + (1:ks)/n^2
amat <- diag(avec)

sigmat1 <- cbind(amat,bvec)
sigmat2 <- cbind(t(bvec),2*n^(-2))
sigmat3 <- rbind(sigmat1,sigmat2)
sigmat <- t(lmat) %*% sigmat3 %*% lmat


# mu vector in Theorem 5@
mubeta <- 1 + 0.5*beta * (beta-1) * diag(sigmat)

#@ Sigma matrix vector in Theorem 5@
sigbeta <- matrix(0,nrow=nrow(sigmat),ncol=ncol(sigmat))
for( i in 1:nrow(sigmat))
{   for( j in 1:ncol(sigmat))
    {sigbeta[i,j] <- beta[i]*beta[j]*sigmat[i,j]}
}

#@ Sum stat in equation 15 @
stat1 <- sum(Vp-1)

#@ QP stat in equation 16 @
QP <- t(Vp-mubeta) %*% solve(sigbeta) %*% (Vp-mubeta)
crit <- qchisq(c(0.01,0.02,0.05,0.10,0.20),df=length(kvec),lower.tail=FALSE)
return(list(Holding.Period=kvec,VRsum=stat1,QPn=QP,ChiSQ.Quantiles_1_2_5_10_20_percent=crit))
}
