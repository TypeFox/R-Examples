init2.k <-
function(reference,response,L,K,t2,m2,previousz,previousclust, previous.alpha, previous.beta,mnr){
x = reference
y = response
maxnr = mnr
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


rdd<-1
if (t2 < 1){
print(paste("exhaustive search is implemented"))
clust<-previousclust
index1<-1:n
cc<-numeric(K-1)
for(k in 1:(K-1)) cc[k]<-length(which(clust==k))
full.components<-which(cc>1)
#big.cluster<-order(cc)[K-1]
m<-length(full.components)
rdd<-0

}

if (t2>0){
m = t2
}





thresh <- -744
if (floor(m)!= m) stop("number of iterations should be positive integer")
if (m < 1) stop("number of iterations should be positive integer")
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
ll<-numeric(m)

# MLE of the library size effects.

for (j in 1:q){
s <- 0
for (l in 1:L[j]){
s <- s+ log(sum(y[,index[j]-1+l]))
}
s <- s/L[j]
for (l in 1:L[j]){
gamma[j,l] <- log(sum(y[,index[j]-1+l])) - s
}}

#number of em runs
ggg<-m2


iter <- 1
for (iter in 1:m){



clust<-previousclust
print(paste("small run: ", iter))
index1<-1:n
if (rdd == 1){
cc<-numeric(K-1)
for(k in 1:(K-1)) cc[k]<-length(which(clust==k))
full.components<-which(cc>1)
#big.cluster<-order(cc)[K-1]
big.cluster<-full.components[floor(runif(1)*length(full.components)+1)]
index1<-which(clust==big.cluster)
}
if(rdd==0){
big.cluster<-full.components[iter]
index1<-which(clust==big.cluster)
}

#this is for random selection of observations within the splitted component
#sam<-sample(length(which(clust==big.cluster)),floor(3+runif(1)*(length(which(clust==big.cluster))-3)))
#sam<-index1[sam]
sam<-index1
#print(paste("random splitting component", big.cluster, ": removing", length(sam), "observations"))
print(paste("splitting component", big.cluster))




z[,1:(K-1)] <- previousz
u.numbers<-runif(length(which(clust==big.cluster)))
z[sam,K]<-z[sam,big.cluster]*u.numbers
z[sam,big.cluster]<-z[sam,big.cluster]*(1-u.numbers)



epsilon <- 1e-10
sl<-length(z[z < epsilon])
bl<-length(z[z > 1-epsilon])
z[z<epsilon]<-rep(epsilon,sl)
z[z>1-epsilon]<-rep(1-epsilon,bl)
z<-z/rowSums(z)


alpha[iter,,1:(K-1)]<-previous.alpha + 0.00*rnorm(length(previous.alpha))
beta[iter,1:(K-1),] <- previous.beta + 0.00*rnorm(length(previous.beta))
for (k in 1:K) {psim[iter,k]  <- sum(z[,k])/n}

for(k in c(big.cluster,K)){
theta<-numeric(q+tau)
# print(k)
nrthreshold <- log(10^(-10));
#maxnr = 10
sc <- nrthreshold + 1
#theta[1:q] <- 2*runif(q)-4# the next elements are the a_{jk}
for (j in 1:q){
sumy<-rowSums(array(y[sam,index[j]:(index[j]-1+L[j])],dim=c(length(sam),1)))
theta[c(j,(q+1):(q+tau))]<-theta[c(j,(q+1):(q+tau))] + as.numeric(glm(sumy~x[sam,],family=poisson,weights = z[sam,k]/previousz[sam,big.cluster])$coefficients)
}
#theta[q+1] <- runif(tau) #beta[iter, k]# the last element is b_{k}
theta[(q+1):(q+tau)]<-theta[(q+1):(q+tau)]/q
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

##################################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################





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




for (gr in 1:ggg){
############################################################################################################################
#E-Step: Mean allocation vectors   #
############################################################################################################################
#print(c("gr = ",gr))
z<-matrix(data = log(psim[iter,]),nrow=n,ncol=K,byrow=T)
for (k in 1:K){
z[,k] <- z[,k] + rowSums(dpois(y,condmeans[[k]],log = T))
}


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
z<-z/rowSums(z)





for(k in 1:K){
# print(k)
nrthreshold <- log(10^(-10));
#maxnr = 10
sc <- nrthreshold + 1
theta[1:q] <- alpha[iter,,k]# the next elements are the a_{jk}
theta[(q+1):(q+tau)] <- beta[iter, k,]# the last element is b_{k}
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
#lll<-mylogLikePoisMix(y, condmeans, psim[iter,])
#print(lll)

}


#bic <- 0
#for(i in 1:n){
#ef <- log(psim[iter,])
#for (j in 1:q) {
#u<-numeric(K);for (t in 1:tau){u<-u + beta[iter,]*x[i,t]}
#for (l in 1:L[j]){
#ef <- ef + dpois(y[i,index[j]-1+l],exp(alpha[iter,j,] + gamma[j,l]  + u),log = T)
#}
#}
#if (max(ef) < thresh) {ef <- c(rep(exp(thresh),K)) } else {ef<-exp(ef)}
#bic <- bic + log(sum(ef)) 
#}


lll<-mylogLikePoisMix(y, condmeans, psim[iter,])


#########################################################################################################################
ll[iter]<-lll$ll
print(ll[iter])
}
print("end of initialization")
ss<- 1:m
mm<-max(ll)
max.index<- ss[ll==mm]
max.index<-max.index[1]
results<-list(psim[max.index,],alpha[max.index,,], beta[max.index,,],ll[max.index])
names(results)<-c("psim","alpha,","beta","ll")
return(results)
}
