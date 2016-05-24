bjmodel <-
function(reference,response,L,m,K,nr,maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z,prev.clust,start.type,prev.alpha,prev.beta){
   
x = reference
y = response
#n.init = initial
n1<-dim(y)[1]
n2<-dim(x)[1]
if (n1 != n2) stop("number of observations does not coincide for x and y")
n = n1
tau<-dim(x)[2]
q<-length(L)
qq<-sum(L)
if (qq != dim(y)[2]) stop("the L vector does not match to the number of columns of y")

index<-numeric(q)
index[1]<-1
if (q > 1){
index[2:q] <- 1+cumsum(L[1:(q-1)])
}
for(j in 1:q){
if (floor(L[j]) != L[j]) stop("L should contain integers > 1")
}

if (K < 1) stop("K should be positive integer")


if (floor(m)!= m) stop("m should be positive integer")

if (m < 1) stop("m should be positive integer")





psim <- matrix(numeric((K*m)),nrow = m,ncol = K,byrow = T);# (m times K) matrix for the weights
z <- matrix(numeric((K*n)),nrow = n,ncol = K,byrow = T);# (n times K) matrix for the allocations
beta <- array(data = 0, dim = c(m,q,tau))# (m times q times tau) matrix for the regression coefficients
alpha <- array(data = NA, dim = c(m,q,K))# (m times q times K) array for the constant terms 
lambda <- numeric(K)# the Lagrange multipliers
theta <-numeric(K + tau + 2)# for the newton raphson iterations
grad <-numeric(K + tau + 2)# the gradient vector
hessian <- array(data = 0, dim = c(K + tau + 2,K + tau + 2))# the hessian matrix
mu <- array(data = NA, dim = c(n,max(L),K))
laa <- numeric(K)
lbb <- numeric(1)
lab <- numeric(K)
bics<-numeric(m)
gamma <- array(data = 0,dim = c(q,max(L)))
delta <- array(data = 0,dim = c(m,q))


# MLE of the library size effects.

for (j in 1:q){
s <- 0
for (l in 1:L[j]){
s <- s+ log(sum(y[,index[j]-1+l]))
}
s <- s/L[j]
for (l in 1:L[j]){
gamma[j,l] <- log(sum(y[,index[j]-1+l])) - s
}
}



iter <- 1

conds = L
n.comp<-K








###########################################################################


# start.type 2 means that we will begin from the previous run
if (start.type == 2) {

print("small-em with splitting initialization scheme")

previousz <- prev.z
previousclust <- prev.clust
previous.alpha<-prev.alpha
previous.beta<-prev.beta
start<-init2.jk.j(reference = x,response = y,L = conds,K = n.comp,tsplit,model = 2,msplit,previousz,previousclust,previous.alpha,previous.beta,mnr = maxnr)



########
} else {

print("2 step small-em with random initialization scheme")


start<-init1.2.jk.j(reference = x,response = y,L = conds,K = n.comp,m1 = m1,m2 = m2,t1 = t1,t2 = t2,model = 2,mnr = maxnr)


}



##########################################################################
print(paste("~~iteration~~","~~~~BIC~~~~","~~loglike diff~~","~~loglike~~"))



psim[iter,]<-start$psim
alpha[iter,,]<-start$alpha
beta[iter,,]<-start$beta
delta[iter,]<-start$delta


# compute the current means
condmeans = vector("list",length = K)
ar<-array(data = NA, dim =c(n,qq))
for (k in 1:K){
i<-0
for (j in 1:q){
#i<-i+1
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + delta[iter,j] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}





theta1<-array(data = NA, dim = c(q,K+2+tau))
for(j in 1:q){
theta1[j,2:(K+1)] <- alpha[iter ,j,]# alpha[j,]'s
theta1[j,K+2] <- delta[iter, j]# delta_{j}
theta1[j,(K+3):(K+2+tau)] <- beta[iter, j,]# beta_{j}
}


#psim[iter,]<-rep(1/K,K)
#for(i in 1:K)psim[iter,k]<-sum(z[,k])/n

d <- q*(K+tau) + K- 1 + sum(L-1)
thresh <- -744
#nrthreshold <-  log(10**(-307))
nrthreshold <-  nr
cl<-10#langrage scale
v<-numeric(K);
emthreshold <-10^(-6)
bics[1] <- -2*start$ll+d*log(n)
criterion <- 10^8

#loglikelist<-numeric(m)

while (criterion >emthreshold&iter<m) {
iter <-iter + 1
############################################################################################################################
#E-Step: Mean allocation vectors   #
############################################################################################################################
z<-matrix(data = log(psim[iter-1,]),nrow=n,ncol=K,byrow=T)
for (k in 1:K){
for (j in 1:q){
for (l in 1:L[j]){
z[,k] <- z[,k] + dpois(y[,index[j]-1+l],condmeans[[k]][,index[j]-1+l],log = T)
}}}


if (K>1){
v1<-which(apply(z,1,max)< thresh)
v3<-1:n
len<-length(v1)
if(len>0){
v2<-apply(array(z[v1,],dim=c(len,K)),1,order)[K,]
ddd<-cbind(v1,v2)
z[v1,]<- 0
z[ddd]<- 1
v3<- -v1 
}
z[v3,]<-exp(z[v3,])}
z<-z/apply(z,1,sum)

epsilon <- 1e-10
sl<-length(z[z < epsilon])
bl<-length(z[z > 1-epsilon])
z[z<epsilon]<-rep(epsilon,sl)
z[z>1-epsilon]<-rep(1-epsilon,bl)
z<-z/apply(z,1,sum)

#q independent maximizations via the newton raphson method



for (j in 1:q){
# initial values
sc <- nrthreshold + 1
theta[1] <- 20*runif(1)# lagrange multiplier 
theta[2:(K+1)] <- alpha[iter - 1,j,]# alpha[j,]'s
theta[K+2] <- delta[iter - 1, j]# delta_{j}
theta[(K+3):(K+2+tau)] <- beta[iter - 1, j,]# beta_{j}
theta1<-theta
#iter1<-0
# Main Loop of N-R iterations
metritis<-1
while (sc > nrthreshold&metritis<maxnr ) {
metritis = metritis + 1
grad <- numeric(K + tau + 2)
#grad[1]<-sum(theta[2:(K+1)])
hessian <- array(data = 0, dim = c(K + tau + 2,K + tau + 2))
# partial derivatives wrt Lagrange multiplier and alpha_{j,}
hessian[1,2:(K+1)]<-1;
for (l in 1:L[j]){
for (k in 1:K){
# computing \mu_{i,l,k} for the current j
mu[,l,k] <- numeric(n)
for (t in 1:tau){
                           mu[,l,k] <- mu[,l,k] + theta[K+2+t]*x[,t];
}
mu[,l,k] <- exp(mu[,l,k] + theta[k+1] + gamma[j,l] + theta[K+2])
#mu[,l,k][mu[,l,k]>5000000]<-50000*runif(length(mu[,l,k][mu[,l,k]>5000000]))
# first derivatives wrt to alpha_{jk}
grad[1+k] <- grad[1+k] + sum(z[,k]*(y[,index[j]-1+l] - mu[,l,k]))  + theta[1]/L[j]
# second derivative wrt to alpha_{jk}
hessian[1+k,1+k] <- hessian[1+k,1+k] - sum(z[,k]*mu[,l,k])
# first derivatives wrt to delta_{j}
grad[K+2] <- grad[K+2] + sum(z[,k]*(y[,index[j]-1+l] - mu[,l,k]))
# second derivatives wrt to delta_{j}
hessian[K+2,K+2] <- hessian[K+2,K+2] - sum(z[,k]*mu[,l,k])
for(t in 1:tau){
# first derivative wrt to beta_{j,t}
grad[K+2+t] <- grad[K+2+t] + sum(z[,k]*(y[,index[j]-1+l] - mu[,l,k])*x[,t])
# second derivative wrt to beta_{j,t}
hessian[K+2+t,K+2+t] <- hessian[K+2+t,K+2+t] - sum(z[,k]*mu[,l,k]*(x[,t]**2))
# partial derivative wrt to alpha_{jk} and beta_{jt}
hessian[1+k,K+2+t] <- hessian[1+k,K+2+t] - sum(z[,k]*mu[,l,k]*x[,t])
# partial derivative wrt to delta_{j} and beta_{jt}
hessian[K+2,K+2+t] <- hessian[K+2,K+2+t] - sum(z[,k]*mu[,l,k]*x[,t])
for (r in seq(t+1,tau,length.out=tau-t)){
# partial derivatives for beta_{jt} and beta_{jr}, r=t+1,...,tau
hessian[K+2+t,K+2+r] <- hessian[K+2+t,K+2+r] - sum(z[,k]*mu[,l,k]*x[,t]*x[,r])
}
}
# partial derivative wrt to alpha_{jk} and delta_{j}
hessian[1+k,K+2] <- hessian[1+k,K+2] - sum(z[,k]*mu[,l,k])
}
}
diag(hessian)<-diag(hessian)/2
hessian<-hessian + t(hessian)
sc <- sum(log(grad[2:(K + tau + 2)]**2))

if (is.nan(sc)==T) 
#{theta<-runif(K+2+tau);sc<-35} else 
{theta<-theta1;sc<-35} else 
if (is.infinite(sc)==T)
#{theta<-runif(K+2+tau);sc<-35} else
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T) theta <- theta1  #theta<-runif(K+2+tau)

#theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
}
alpha[iter,j,] <- theta[2:(K+1)]# alpha[j,]'s
delta[iter, j] <- theta[K+2]# delta_{j}
beta[iter, j,] <- theta[(K+3):(K+2+tau)]# beta_{j}

}
#############################################
#Maximizing according to the weights #
#############################################
for (k in 1:K) {psim[iter,k]  <- sum(z[,k])/n}

# compute the current means
condmeans = vector("list",length = K)
ar<-array(data = NA, dim =c(n,qq))
for (k in 1:K){
i<-0
for (j in 1:q){
#i<-i+1
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + delta[iter,j] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}

lll<-mylogLikePoisMix(y, condmeans, psim[iter,])
#loglikelist[iter]<-lll[[1]]
#plot(loglikelist[2:iter],type="l",main="loglikelihood",xlab="iteration")
bic<- -2*lll[[1]]+d*log(n)
bics[iter] <- bic
criterion <-  (abs(bics[iter - 1] - bics[iter]))  #/bics[iter - 1]
print(c(iter,bic,criterion,lll[[1]]))

}




#print(iter)
clust <- numeric(n)
for(i in 1:n) clust[i] <- order(z[i,])[K]

nz<-z
ind<-1:K

for (i in 1:n){
index<-ind[nz[i,]< exp(thresh)]
nz[i,index]<- rep(exp(thresh),length(index))
nz[i,]<- nz[i,]/sum(nz[i,])
}



iclbic <- bic - 2*sum(nz*log(nz))
print(c("ICL = ", iclbic))



for(j in 1:q){
for(k in 1:K){
alpha[1:iter,j,k]<-alpha[1:iter,j,k] + delta[1:iter,j]
}}
alpha<-array(data=alpha[1:iter,,],dim=c(iter,q,K))
beta <- array(data=beta[1:iter,,],dim=c(iter,q,tau))
psim<-psim[1:iter,]


results<-list(alpha,beta,gamma,psim,clust,z,bic,iclbic,lll[[1]])
names(results)<-c("alpha","beta","gamma","psim","clust","z","bic","icl","ll")
return(results)


}
