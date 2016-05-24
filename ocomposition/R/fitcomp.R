fitcomp <- function(data.v, data.x, n.formula, v.formula, l.bound=1, n.sample=100, burn=0, thin=1, init=NULL){

cat("ocomposition package version 1.1", '\n')
cat(" ",'\n')

v=data.matrix(data.v)
v[is.na(v)]=0
if (sum(round(apply(v,1,sum),5))/nrow(v)!=1) stop("Components don't add up to 1. Check your data or rescale")
v=t(apply(v,1,sort,decreasing=TRUE))
t = which(apply(v, 2, sd) == 0)
if (length(t) > 0) {v = v[, -t]}

n=apply(v,1,function(x) sum(x > 0))
x=data.frame(data.x, n=n)
y=t(apply(v,1,v.y)) # map compsitions into R space
rownames(y) = names(n) = rownames(data.v)

X=model.matrix(n.formula, data=x)
G=model.matrix(v.formula, data=x)

d = intersect(intersect(rownames(X), rownames(G)), rownames(y))
y = y[d, ]
X = X[d, ]
G = G[d, ]
n = n[d]

if (any(is.na(X))) stop("Missing values in the covariates are not allowed.")

if (sd(n, na.rm=TRUE)==0) stop("The number of components does not vary. n cannot be used as a covariate. Change v.formula accordingly.")

                     

N=length(n); 
D=max(n)-1; # dimension of the model
K=ncol(X);
P=ncol(G);

#=====================================================
#  COUNT MODEL PARAMETERS 
#====================================================

# L-LIKELIHOOD
n.lik=function(b){
sum(log(dtnegbin(n,mu=exp(X%*%b[1:K]),dispersion=exp(b[K+1]),l.bound=l.bound)))
}


if (sd(n) > 0){
# MLE ESTIMATES FOR TRUNCATED N BIN MODEL
n.mle=optim(rep(0,K+1), n.lik, method="L-BFGS-B",hessian=TRUE, control=list(fnscale=-1),lower=c(rep(-Inf,K),-5),upper=c(rep(Inf,K),5))
H=-0.2*solve(n.mle$hessian/K)
}

if (sd(n) ==0) {
	n.mle=list(par=rep(0, K))
	cat("Component counts do not vary across the units. Count model omitted.", '\n')
}



# LOG POSTERIOR
pois.post=function(b){
sum=sum(log(dtnegbin(n, mu = exp(X%*%b$fix[1:K]), dispersion=exp(b$fix[K+1]), l.bound=l.bound))) + dnorm(b$fix[K+1], 0, sqrt(10), log=TRUE) 
sum
}


# SAMPLING FROM LOG-POSTERIOR
update.b=function(b){
cand = b
cand$fix=c(rmnorm(1,mean=b$fix,varcov=H))
prob=min(0,pois.post(cand)-pois.post(b))
if (runif(1) < exp(prob)){b$fix = cand$fix}
b
}


#========================================================# 		
# COMPOSITIONAL PARAMETERS
#========================================================

# Z, gamma, mu, rho

#========================================================
#---     Sweep operator for conditional MVN sampling
#========================================================

cov.f=function(x,Sigma){
if (x[D+1]==D){return(x[-(D+1)])}
if (x[D+1]<D){
l=1:x[D+1]
t=(length(l)+1):D
S=Sigma[t,t]-Sigma[t,l]%*%solve(Sigma[l,l])%*%Sigma[l,t]
m=c(x[l],x[t]%*%chol(S))
return(m)
}
}

mvm=function(x,Sigma,g){
if (x[P+D+1] == D){return(rep(0,D))}
if (x[P+D+1] < D){
l=1:x[P+D+1]
t=(length(l)+1):D
m=x[1:P]%*%g
M=m[t]+Sigma[t,l]%*%solve(Sigma[l,l])%*%(x[(P+1):(P+length(l))]-m[l])
return(c(rep(0,length(l)),M))
}
}


# SAMPLE LATENT PARAMETER Z
update.Z=function(y,n,g,Sigma,tau){
G.n=cbind(G,y,n-1)
Z=rmnorm(N,0,varcov=diag(D))/sqrt(tau)
X.n=cbind(Z,n-1)
M=t(apply(G.n,1,mvm,Sigma=Sigma,g=g))
Z=t(apply(X.n,1,cov.f,Sigma=Sigma))+M
Z[!is.na(y)]=y[!is.na(y)]
return(Z)
}


# SAMPLE REGRESSION COEFFICIENTS FOR Y MODEL
update.g=function(Z,Sigma,mu,rho,tau){
sum.xx=kronecker(solve(Sigma),t(G*tau)%*%G)
sum.xy=kronecker(solve(Sigma),t(G*tau))%*%c(Z)
G.prior=rep(rho,each=D)*diag(P*D)
g.prior=rep(mu,each=D)
V=solve(sum.xx+solve(G.prior))
M=V%*%(sum.xy+solve(G.prior)%*%g.prior)
g=rmnorm(1,mean=M,varcov=V)
return(matrix(g,P,D))
}


# HYPERPARAMETERS FOR COEFFICIENTS
update.mu=function(g,rho){
mean=apply(g,1,mean)
mu=rnorm(P,mean=mean,sd=sqrt(rho/D))
return(mu)
}

update.rho=function(g,mu,al,be){
W=apply((g-mu)^2,1,sum)
return(1/rgamma(P,shape=D/2+al,rate=W/2+be))
}

# Inverse-Wishart samples
 riwish = function(v0, S0) {
 X <- rmnorm(n=round(v0), varcov=solve(S0))
 solve(t(X)%*%X)
 }

# SAMPLE COVARIANCE MATRIX
update.Sigma=function(Z,g,psi,tau){
e=sqrt(tau)*(Z-G%*%g)	
v_0=D+1+exp(psi[1]);  B=diag(D)*exp(sum(psi)) + t(e)%*%(e)
return(riwish(v_0 + N, B))
}

# SAMPLE HYPERPARAMETERS OF THE COV MATRIX
wpost=function(psi,Sigma,sd){
v_0=D+1+exp(psi[1])
IW=diag(D)*(exp(sum(psi)))
lpost=lndIWishart(v_0, IW,Sigma) + sum(dnorm(c(psi),0,sd,log=TRUE))
return(lpost)
}

update.psi=function(step,psi,Sigma,sd){
cand=c(psi+rnorm(2,0,step))
prob=min(0,wpost(cand,Sigma,sd)-wpost(psi,Sigma,sd))
if(runif(1)<exp(prob)){psi=cand}
return(psi)
}

# SAMPLE LATENT MIXING PARAMETER TAU
m.dist=function(x,Sigma){t(x)%*%solve(Sigma)%*%x}
update.tau=function(Z,nu,g,Sigma){
rgamma(N,shape=(nu+D)/2,rate=apply(Z-G%*%g,1,m.dist,Sigma=Sigma)/2+nu/2)
}


# SAMPLE DEGREES OF FREEDOM PARAMETER NU

lpost.nu=function(nu,tau){sum(dgamma(tau, shape=nu/2, rate=nu/2,log=TRUE))-2*log(nu)}
update.nu=function(nu,tau,step){
cand=1+rgamma(1,shape=nu*step,rate=step)
prob=min(0,lpost.nu(cand,tau)-lpost.nu(nu,tau)-dgamma(cand-1,shape=nu*step,rate=step,log=TRUE)+dgamma(nu,shape=nu*step,rate=step,log=TRUE))
if (runif(1) < exp(prob)){nu = cand}
return(nu)
}


#### STARTING VALUES

if(is.null(init)){
out=list(
Z=NULL,
b=list(fix=n.mle$par, sigma2=1), 
g=rmnorm(P,0,.5*diag(D)),
mu=rep(0,P),
rho=rep(1,P),
Sigma=diag(D),
psi=c(0,0),
tau=rep(1,N),
nu=10
)
}

if(!is.null(init)) out=init

samples=list(
g=array(NA, c(P,D,n.sample)),
b=matrix(NA, nrow=n.sample, ncol=K+1),
rand.sigma2 = rep(NA, n.sample),
mu=matrix(NA, nrow=n.sample, ncol=P),
rho=matrix(NA, nrow=n.sample, ncol=P),
psi=matrix(NA,nrow=n.sample, ncol=2),
Sigma=array(NA,c(D,D,n.sample)),
tau=matrix(NA,n.sample,N),
nu=rep(NA,n.sample),
D=D,
l.bound=l.bound,
n.formula=n.formula,
v.formula=v.formula
)

colnames(samples$b)=c(colnames(X),"ln(w)")
colnames(samples$mu)=colnames(samples$rho)=colnames(G)
colnames(samples$g)=c(paste("y",1:D,sep="-"))
rownames(samples$g)=colnames(G)

gibbs=function(out){
  
out$Z=update.Z(y,n,out$g,out$Sigma,out$tau)
out$tau=update.tau(out$Z,out$nu,out$g,out$Sigma)
out$nu=update.nu(out$nu,out$tau,step=1.5)
out$g=update.g(out$Z,out$Sigma,out$mu,out$rho,out$tau)
if (sd(n) > 0) out$b=update.b(out$b)
out$mu=update.mu(out$g,out$rho)
out$rho=update.rho(out$g,out$mu,al=1e-6,be=1e-6)
out$Sigma=update.Sigma(out$Z,out$g,out$psi,out$tau)
out$psi=update.psi(.4,out$psi,out$Sigma,sd=1)
return(out)
}


n.iter = burn + n.sample*thin
pars = length(unlist(out)) - (D-1)*D/2

cat("Number of iterations:", n.iter, '\n')
cat("Number of parameters:", pars,'\n')
cat("Number of units:     ", N, '\n')
main.time <- proc.time()
out=gibbs(out) 
cat("Estimated time:      ", round(n.iter*((proc.time()-main.time)/60)[1],2) + 0.0014*n.iter, "mins", '\n')

cat("                    ", '\n')
cat("Progress: ")

t=1
while(t<=burn) { ### burn-in chain
out=gibbs(out) 
if (t%%(n.iter/10)==0)  cat(100*t/n.iter, "%  ", sep="")
t=t+1
}


t=1
while(t<=n.sample){
for (i in 1:thin){
out=gibbs(out)
samples$Sigma[,,t]=out$Sigma
samples$g[,,t]=out$g
if(sd(n) > 0) {
  samples$b[t,]=unlist(out$b$fix)
}
  samples$mu[t,]=out$mu
samples$rho[t,]=out$rho
samples$psi[t,]=out$psi
samples$tau[t,]=out$tau
samples$nu[t]=out$nu
}
if ((t+burn)%%(n.iter/10)==0){
cat(100*(t+burn)/n.iter, "%  ", sep="")
}
t=t+1
}

cat("",'\n')

attr(samples, "date") <- date()
attr(samples$b, "ext.name") <- "Component count model: beta parameters"
attr(samples$mu, "ext.name") <- "Component size model hyper-parameter: mu"
attr(samples$rho, "ext.name") <- "Shrinkage hyper-parameters: rho"

class(samples$b) = class(samples$g) = class(samples$mu) = class(samples$rho) = "out"

class(samples) = "composition"

cat('\n')
cat("Time elapsed:", round(((proc.time()-main.time)/60)[1],2), "mins",'\n')

samples

}