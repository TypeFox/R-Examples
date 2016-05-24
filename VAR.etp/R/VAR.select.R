VAR.select <-
function(x,type="const",ic="aic",pmax)
{
K <- ncol(x)
icmat <- matrix(NA,nrow=pmax+1)
m=t(t(x)-colMeans(x))

for(i in 0:pmax)
{
T <- nrow(x)-i
if(i==0) sigu=(t(m) %*% m)/nrow(x)
if(i > 0) sigu <- VAR.est(x,p=i,type)$sigu *( (T-K*i-1)/(T))
LL <- log(det(sigu)); n <- i*K^2
#if(type=="none") n <- i*K^2
#if(type=="const") n <- i*K^2 + K
#if(type=="const+trend") n <- i*K^2 + 2*K

if(ic =="aic") icmat[i+1] <- LL + 2*n/(T)
if(ic == "sc") icmat[i+1] <- LL + log(T)*n/(T)
if(ic == "hq") icmat[i+1] <- LL + 2*log(log(T))*n/(T)

}
ps <- which.min(icmat)-1
colnames(icmat) <- ic; rownames(icmat) <- 0:pmax
return(list(IC=icmat,p=ps))
}
