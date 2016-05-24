bjkmodel <-
function(reference,response,L,m,K,nr,maxnr,m1,m2,t1,t2, msplit,tsplit,prev.z,prev.clust,start.type,prev.alpha,prev.beta){
x = reference
y = response
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




psim <- matrix(numeric((K*m)),nrow = m,ncol = K,byrow = T);# (m times K) matrix for the weights
z <- matrix(numeric((K*n)),nrow = n,ncol = K,byrow = T);# (n times K) matrix for the allocations
beta <- array(data = 0, dim = c(m,q,K,tau))# (m times q times tau) matrix for the regression coefficients
alpha <- array(data = NA, dim = c(m,q,K))# (m times q times K) array for the constant terms 
lambda <- numeric(K)# the Lagrange multipliers
theta <-numeric(tau + 1)# for the newton raphson iterations
grad <-numeric(tau + 1)# the gradient vector
hessian <- array(data = 0, dim = c(tau + 1,tau + 1))# the hessian matrix
mu <- array(data = NA, dim = c(n,max(L)))
laa <- numeric(K)
lbb <- numeric(1)
lab <- numeric(K)
bics<-numeric(m)
gamma <- array(data = 0,dim = c(q,max(L)))


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



# start.type 2 means that we will begin from the previous run
if (start.type == 2) {

print("small-em with splitting initialization scheme")

previousz <- prev.z
previousclust <- prev.clust
previous.alpha<-prev.alpha
previous.beta<-prev.beta
start<-init2.jk.j(reference = x,response = y,L = conds,K = n.comp,tsplit,model = 1,msplit,previousz,previousclust,previous.alpha,previous.beta,mnr = maxnr)

} else {

print("2 step small-em with random initialization scheme")

start<-init1.2.jk.j(reference = x,response = y,L = conds,K = n.comp,m1=m1,m2=m2,t1=t1,t2=t2,model = 1,mnr = maxnr)


}



psim[iter,]<-start$psim
alpha[iter,,]<-start$alpha
beta[iter,,,]<-start$beta

print(paste("~~iteration~~","~~~~BIC~~~~","~~loglike diff~~","~~loglike~~"))


# compute the current means
condmeans = vector("list",length = K)
ar<-array(data = NA, dim =c(n,qq))
for (k in 1:K){
i<-0
for (j in 1:q){
#i<-i+1
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,k,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}





d <- q*K*(1+tau) + K- 1 + sum(L-1)
thresh <- -744
#nrthreshold <-  log(10**(-307))
#nrthreshold <-  log(10**(-100))/10
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
z[,k] <- z[,k] + rowSums(dpois(y,condmeans[[k]],log = T))
}
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
z[v3,]<-exp(z[v3,])
z<-z/rowSums(z)

epsilon <- 1e-10
sl<-length(z[z < epsilon])
bl<-length(z[z > 1-epsilon])
z[z<epsilon]<-rep(epsilon,sl)
z[z>1-epsilon]<-rep(1-epsilon,bl)

}
z<-z/rowSums(z)




#q independent maximizations via the newton raphson method


for (k in 1:K){
for (j in 1:q){
# initial values
sc <- nrthreshold + 1
theta[1] <- alpha[iter - 1,j,k]# alpha[j,]'s
theta[2:(tau+1)] <- beta[iter - 1, j,k,]# beta_{j}
# Main Loop of N-R iterations
metritis<-1
while (sc > nrthreshold&metritis<maxnr ) {
metritis = metritis + 1
grad <- numeric(tau + 1)
hessian <- array(data = 0, dim = c(tau + 1,tau + 1))
for (l in 1:L[j]){
# computing \mu_{i,l,k} for the current j
mu[,l] <- numeric(n)
for (t in 1:tau){
                           mu[,l] <- mu[,l] + theta[1+t]*x[,t];
}
mu[,l] <- exp(mu[,l] + theta[1] + gamma[j,l])
# first derivatives wrt to alpha_{jk}
grad[1] <- grad[1] + sum(z[,k]*(y[,index[j]-1+l] - mu[,l]))
# second derivative wrt to alpha_{jk}
hessian[1,1] <- hessian[1,1] - sum(z[,k]*mu[,l])
for(t in 1:tau){
# first derivative wrt to beta_{j,t}
grad[1+t] <- grad[1+t] + sum(z[,k]*(y[,index[j]-1+l] - mu[,l])*x[,t])
# second derivative wrt to beta_{j,t}
hessian[1+t,1+t] <- hessian[1+t,1+t] - sum(z[,k]*mu[,l]*(x[,t]**2))
# partial derivative wrt to alpha_{jk} and beta_{jt}
hessian[1,1+t] <- hessian[1,1+t] - sum(z[,k]*mu[,l]*x[,t])
for (r in seq(t+1,tau,length.out=tau-t)){
# partial derivatives for beta_{jt} and beta_{jr}, r=t+1,...,tau
hessian[1+t,1+r] <- hessian[1+t,1+r] - sum(z[,k]*mu[,l]*x[,t]*x[,r])
}
}
}
diag(hessian)<-diag(hessian)/2
hessian<-hessian + t(hessian)
sc <- sum(log(grad**2))
#if (is.infinite(sc)==T) {theta<-runif(1+tau);print("oops")} else
#theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
if (is.nan(sc)==T) 
#{theta<-runif(1+tau);sc<-35} else 
{theta[1]<-alpha[iter-1,j,k];theta[2:(tau+1)]<-beta[iter-1,j,k,];sc<-35} else 
if (is.infinite(sc)==T)
#{theta<-runif(1+tau);sc<-35} else
{theta[1]<-alpha[iter-1,j,k];theta[2:(tau+1)]<-beta[iter-1,j,k,];sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T) {theta[1]<-alpha[iter-1,j,k];theta[2:(tau+1)]<-beta[iter-1,j,k,]}#  theta<-runif(1+tau) #   #

}
#print(metritis)
alpha[iter,j,k] <- theta[1]# alpha[j,]'s
beta[iter, j,k,] <- theta[2:(1+tau)]# delta_{j}

}}
#############################################
#Maximizing according to the weights #
#############################################
for (k in 1:K) {psim[iter,k]  <- sum(z[,k])/n}

condmeans = vector("list",length = K)
ar<-array(data = NA, dim =c(n,qq))
for (k in 1:K){
i<-0
for (j in 1:q){
#i<-i+1
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,k,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}

lll<-mylogLikePoisMix(y, condmeans, psim[iter,])
#loglikelist[iter]<-lll[[1]]
#plot(loglikelist[2:iter],type="l",main="loglikelihood",xlab="iteration")

bic<--2*lll[[1]]+d*log(n)
bics[iter] <- bic
criterion <-  0.5*abs(bics[iter - 1] - bics[iter]) #/bics[iter - 1]
print(c(iter,bic,criterion,lll[[1]]))

}
#print(iter)

clust <- numeric(n)
for(i in 1:n) clust[i] <- order(z[i,])[K]


# computing icl after smoothing the z's

nz<-z
ind<-1:K

for (i in 1:n){
index<-ind[nz[i,]< exp(thresh)]
nz[i,index]<- rep(exp(thresh),length(index))
nz[i,]<- nz[i,]/sum(nz[i,])
}



iclbic <- bic - 2*sum(nz*log(nz))
print(c("ICL = ", iclbic))



alpha<-array(data=alpha[1:iter,,],dim = c(iter,q,K))
beta <- array(beta[1:iter,,,],dim = c(iter,q,K,tau))
psim<-psim[1:iter,]

#return(alpha)
#return(beta)
#return(gamma)
#return(psim)
#return(clust)

results<-list(alpha,beta,gamma,psim,clust,z,bic,iclbic,lll[[1]])
names(results)<-c("alpha","beta","gamma","psim","clust","z","bic","icl","ll")
return(results)


}
