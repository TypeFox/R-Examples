Bayesian.SGoF <-
function(u,alpha=0.05,gamma=0.05,P0=0.5,a0=1,b0=1){

bayesian.sgof<-function(u,alpha=0.05,gamma=0.05,P0=0.5,a0=1,b0=1){


n=length(u)


ss=length(u[u<=gamma])
P1<-1-P0


s=seq(0,n)


p=gamma;rho=seq(0.001,0.999,by=0.001)
PHoxx<-matrix(nrow=n+1,ncol=length(rho))
BBxx<-matrix(nrow=n+1,ncol=length(rho))
ax=(1-rho)*p/rho;bx=(1-rho)*(1-p)/rho
for(j in 1:(n+1)){
BBxx[j,]=exp(log(P1)-log(P0)-s[j]*log(gamma)+(s[j]-n)*log(1-gamma)+lbeta(ax+s[j],n-s[j]+bx)-lbeta(ax,bx))
}
for(j in 1:(n+1)){
PHoxx[j,]=1/(1+((1-P0)/P0)*BBxx[j,])
}



ind<-sapply(1:(n+1), function(j) ind<-which.min(PHoxx[j,]))
lowb<-sapply(1:(n+1), function(j) lowb<-PHoxx[j,ind[j]])

if(sum(lowb > 0.05)==0){salpha=0}else{salpha = max(which(lowb > 0.05))}


Bayesian.SGoF = max(min(floor((ss >= salpha) *( n * (qbeta(alpha,a0 + ss, b0 + n - ss) - gamma) + 1)), sum(as.integer(n *ecdf(u)(u)) <= (ss >= salpha) *( n * (qbeta(alpha,a0 + ss, b0 + n - ss) - gamma) + 1))), 0)



theta0=gamma

BB=exp(log(P1)-log(P0)-ss*log(theta0)+(ss-n)*log(1-theta0)+lbeta(a0+ss,n-ss+b0)-lbeta(a0,b0))
PHo=round(1/(1+BB),5)




su<-sort(u)
jj<-which(u==1)
if(length(jj)!=0) pi0<-1 else pi0<-min((-1/n)*sum(log(1-u)),1)

if(Bayesian.SGoF==0){FDR_BS<-0}else{FDR_BS<-round((pi0*su[Bayesian.SGoF])/(ecdf(u)(su[Bayesian.SGoF])),4)}


return(c(list(Rejections=Bayesian.SGoF,FDR=min(FDR_BS,1),Posterior=PHo,s=ss,s.alpha=salpha)))
}





if(missing(u)){stop("data argument is required")}
u<-as.vector(u)
res<-bayesian.sgof(u,alpha,gamma,P0,a0,b0)
res$data<-sort(u)
res$alpha<-alpha
res$gamma<-gamma
res$P0<-P0
res$a0<-a0
res$b0<-b0
res$call<-match.call()
class(res)<-"Bayesian.SGoF"
return(res)
}
