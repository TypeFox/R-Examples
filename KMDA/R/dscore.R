dscore <-
function(x,y,lower,upper,m){  
if (ncol(x)==1){x=t(x)} ##If the set contains only one metabolite, then transpose it to a row vector.
n=ncol(x)
fit=glm(y~1,family=binomial)  ## fit the null model
betahat=fit$coef
mu_0=invlogit(betahat) 
D_0=diag(mu_0*(1-mu_0),n)
X=rep(1,n)
P_0=D_0-D_0%*%X%*%solve(t(X)%*%D_0%*%X)%*%t(X)%*%D_0
rho=seq(lower,upper,length=m)
S = rep(0,m)
for(l in 1:m){
k=matrix(NA,n,n)  
	for (i in 1:n){
		for(j in 1:n){
			k[i,j]=dkernel(x[,i],x[,j],rho[l])
	}
}
Q=t(y-mu_0)%*%k%*%(y-mu_0)
mu= tr(P_0%*% k)  
sigma=sqrt(2*tr(P_0%*% k%*%P_0%*% k))
S[l]=(Q-mu)/sigma
}
M=max(S)
W=0
for(i in 1:(m-1)){
	W=W+abs(S[i+1]-S[i])
}
pvalue=pnorm(-M)+W*exp(-M^2/2)/sqrt(8*pi) ## (12) in Liu et al. (2008)
return(pvalue)
}
