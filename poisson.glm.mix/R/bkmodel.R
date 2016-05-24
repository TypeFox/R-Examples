bkmodel <-
function(reference,response,L,m,K,nr,maxnr,t2,m2,prev.z,prev.clust,start.type,prev.alpha,prev.beta){
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

#if (tau>1)stop("x shoul be one variable")

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
beta <- array(data = NA, dim = c(m,K,tau))# (m times K) matrix for the regression coefficients
alpha <- array(data = NA, dim = c(m,q,K))# (m times q times K) array for the constant terms 
lambda <- numeric(K)# the Lagrange multipliers
theta <-numeric(q+tau)# for the newton raphson iterations
grad <-numeric(q+tau)# the gradient vector
hessian <- array(data = 0, dim = c(q+tau,q+tau))# the hessian matrix
mu <- array(data = NA, dim = c(n,q,max(L)))
laa <- numeric(q)
lbb <- numeric(tau)
lab <- numeric(q*tau)
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

if (start.type==1){
print("small-em with random initialization scheme")
start<-init1.k(reference = x,response = y,L = conds,K = n.comp,t2,m2,mnr = maxnr)
}else{
print("small-em with splitting initialization scheme")

previousz <- prev.z
previousclust <- prev.clust
previous.alpha<-prev.alpha
previous.beta<-prev.beta

start<-init2.k(reference=x,response=y,L=conds,K=n.comp,t2,m2,previousz,previousclust, previous.alpha, previous.beta,mnr = maxnr)

}
psim[iter,]<-start$psim
alpha[iter,,]<-start$alpha
beta[iter,,]<-start$beta


print(paste("~~iteration~~","~~~~BIC~~~~","~~loglike diff~~","~~loglike~~"))




# compute the current means
condmeans = vector("list",length = K)
ar<-array(data = NA, dim =c(n,qq))
for (k in 1:K){
i<-0
for (j in 1:q){
#i<-i+1
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + gamma[j,l] + colSums(beta[iter,k,]*t(x)))
}}
condmeans[[k]]<-ar
}





d <- (q+tau)*K + K - 1 + sum(L-1)
thresh <- -744
#nrthreshold <-  log(10**(-307))
#nrthreshold <-  log(10**(-100))/10
nrthreshold <-  nr
cl<-10#langrage scale
v<-numeric(K);
emthreshold <-10^(-6)
bics[1] <- -2*start$ll+d*log(n)
criterion <- 10^10


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
z[v3,]<-exp(z[v3,])}
z<-z/rowSums(z)


epsilon <- 1e-10
sl<-length(z[z < epsilon])
bl<-length(z[z > 1-epsilon])
z[z<epsilon]<-rep(epsilon,sl)
z[z>1-epsilon]<-rep(1-epsilon,bl)
z<-z/rowSums(z)










#q independent maximizations via the newton raphson method



for(k in 1:K){
sc <- nrthreshold + 1
theta[1:q] <- alpha[iter-1,,k]# the next elements are the a_{jk}
theta[(q+1):(q+tau)] <- beta[iter-1, k,]# the last element is b_{k}
metritis <- 1
theta1<-theta
while (sc > nrthreshold&metritis<maxnr ) {
hessian<-array(data=0,dim=c(q+tau,q+tau))
grad<- numeric(q+tau)
metritis <- metritis + 1
uo <- 0
mudot <- numeric(n)
for (j in 1:q){ 
for (l in 1:L[j]){mu[,j,l] = exp(theta[j] + gamma[j,l] + colSums(theta[(q+1):(q+tau)]*t(x)))# store the means
}
#################################################################
#grad[j] <- 0;
laa[j]<-0;
lab[j]<-0;
for (l in L[j]){
grad[j] <- grad[j] + sum(z[,k]*(y[,index[j]-1+l]-mu[,j,l])); 
for (t in 1:tau) {
grad[q+t] <- grad[q+t] + sum(z[,k]*(y[,index[j]-1+l]-mu[,j,l])*x[,t])
hessian[j,q+t] <- hessian[j,q+t] - sum(z[,k]*mu[,j,l]*x[,t])# second partial derivatives wrt to a_{jk} and b_{kt}
for (r in seq(t+1,tau,length.out=tau-t)){
# partial derivatives for beta_{jt} and beta_{jr}, r=t+1,...,tau
hessian[q+t,q+r] <- hessian[q+t,q+r] - sum(z[,k]*mu[,j,l]*x[,t]*x[,r])
}

}
mudot <- mudot + mu[,j,l];
laa[j] <- laa[j] - sum(z[,k]*mu[,j,l]);# second partial derivatives wrt to a_{jk}
}
#####################################################################
#computing the partial derivatives for component k       #
#####################################################################
}   
for (t in 1:tau) hessian[q+t,q+t] <- -sum(z[,k]*mudot*(x[,t]**2))######### second partial derivatives wrt to b_{kt} 
#################################################
#constructing the hessian matrix#
#################################################
for (j in 1:q){hessian[j,j] <- laa[j]}
diag(hessian)<-diag(hessian)/2
hessian<-hessian + t(hessian)

sc<-sum(log(grad**2))
if (is.nan(sc)==T) 
{theta<-theta1;sc<-35} else  
if (is.infinite(sc)==T)
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
if (is.na(max(theta))==T)theta<-theta1
}
############################################
#updated parameters for component k #
############################################
alpha[iter,,k] <- theta[1:q]
beta[iter, k,] <- theta[(q+1):(q+tau)]
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
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + gamma[j,l] + colSums(beta[iter,k,]*t(x)))
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

alpha<-array(alpha[1:iter,,],dim=c(iter,q,K))
beta <- array(beta[1:iter,,],dim=c(iter,K,tau))
psim<-psim[1:iter,]




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



results<-list(alpha,beta,gamma,psim,clust,z,bic,iclbic,lll[[1]])
names(results)<-c("alpha","beta","gamma","psim","clust","z","bic","icl","ll")
return(results)


}
