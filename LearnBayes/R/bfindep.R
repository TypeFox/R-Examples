bfindep=function(y,K,m)
{
# compute Bayes factor against independence
# using Albert and Gupta independence priors
# ymat - I x J matrix
# K - Dirichlet precision parameter
# m - number of iterations

rdirichlet=function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

ldirichlet=function(alpha)
{
# log dirichlet function 
# for multiple values stored in matrix alpha
return(rowSums(lgamma(alpha))-lgamma(rowSums(alpha)))
}

yc=colSums(y); yr=rowSums(y); n=sum(yc)
d=dim(y); I=d[1]; J=d[2]

etaA=rdirichlet(m,yr+1)
etaB=rdirichlet(m,yc+1)

Keta=c(); KetaY=c()
for (i in 1:I)
{
for (j in 1:J)
{
Keta=cbind(Keta,K*etaA[,i]*etaB[,j])
KetaY=cbind(KetaY,K*etaA[,i]*etaB[,j]+y[i,j])
}}

logint=ldirichlet(KetaY)-ldirichlet(Keta)
for (i in 1:I) logint=logint-yr[i]*log(etaA[,i])
for (j in 1:J) logint=logint-yc[j]*log(etaB[,j])

int=exp(logint)
  
return(list(bf=mean(int),nse=sd(int)/sqrt(m)))
}
