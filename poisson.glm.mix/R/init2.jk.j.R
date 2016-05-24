init2.jk.j <-
function(reference,response,L,K,tsplit,model,msplit,previousz,previousclust, previous.alpha, previous.beta,mnr){
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
start.type = 1
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
if (tsplit < 1){
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

if (tsplit>0){
m = tsplit
}
thresh <- -744
#if (floor(m)!= m) stop("number of iterations should be positive integer")

if (m < 1) stop("number of iterations should be positive integer")

if (model == 1){beta <- array(data = 0, dim = c(m,q,K,tau))}
if (model == 2){beta <- array(data = 0, dim = c(m,q,tau))}

psim <- matrix(numeric((K*m)),nrow = m,ncol = K,byrow = T);# (m times K) matrix for the weights
z <- matrix(numeric((K*n)),nrow = n,ncol = K,byrow = T);# (n times K) matrix for the allocations
alpha <- array(data = NA, dim = c(m,q,K))# (m times q times K) array for the constant terms 
lambda <- numeric(K)# the Lagrange multipliers
theta <-numeric(tau + 1)# for the newton raphson iterations
grad <-numeric(tau + 1)# the gradient vector
hessian <- array(data = 0, dim = c(tau + 1,tau + 1))# the hessian matrix
mu <- array(data = NA, dim = c(n,max(L)))
if (model == 2){ mu <- array(data = NA, dim = c(n,max(L),K)) }
laa <- numeric(K)
lbb <- numeric(1)
lab <- numeric(K)
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




iter <- 1



if (model == 1){
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
if (1<2){

#####################old scheme
#this is for random selection of observations within the splitted component
#sam<-sample(length(which(clust==big.cluster)),floor(3+runif(1)*(length(which(clust==big.cluster))-3)))
#sam<-index1[sam]
sam<-index1
print(paste("splitting component", big.cluster))

} else{
#this is for kmeans splitting
print(paste("cluster=", big.cluster))
cl<-kmeans(log(rowSums( y[clust==big.cluster,]) +0.5 )/rowSums(array(x[clust==big.cluster,]+0.5,dim=c(length(x[clust==big.cluster,]),tau))),2)
len1<-length(index1[cl$cluster==1])
len2<-length(index1[cl$cluster==2])
print(paste("lengths=", len1,len2))
if (len1*len2==0){
sam<-sample(length(which(clust==big.cluster)),floor(1+length(which(clust==big.cluster))/(1+10*runif(1))))
sam<-index1[sam]
} else{

previousz[index1[cl$cluster==2],big.cluster]<-1
sam<- index1[cl$cluster==1]
}
print(paste("kmeans-splitting component", big.cluster, " n = ", length(sam)))

}

z[,1:(K-1)] <- previousz
u.numbers<-runif(length(which(clust==big.cluster)))
#u.numbers<-rbeta(length(which(clust==big.cluster)),0.1,0.1)
z[sam,K]<-z[sam,big.cluster]*u.numbers
z[sam,big.cluster]<-z[sam,big.cluster]*(1-u.numbers)



epsilon <- 1e-10
sl<-length(z[z < epsilon])
bl<-length(z[z > 1-epsilon])
z[z<epsilon]<-rep(epsilon,sl)
z[z>1-epsilon]<-rep(1-epsilon,bl)
z<-z/rowSums(z)


#######################################################################################################################
#######################################################################################################################
###############################initialization of values################################################################
#######################################################################################################################
#######################################################################################################################

alpha[iter,,1:(K-1)]<-previous.alpha + 0.00*rnorm(length(previous.alpha))
beta[iter,,1:(K-1),] <- previous.beta + 0.00*rnorm(length(previous.beta))
for (k in 1:K) {psim[iter,k]  <- sum(z[,k])/n}


for(k in c(big.cluster,K)){
for(j in 1:q){
#st[1]<-alpha[iter,j,k]
#st[2:(tau+1)]<-beta[iter,j,k,]

nrthreshold <- log(10^(-10));
#maxnr = 10
sc <- nrthreshold + 1
#theta<- runif(1+tau)
#sumy<-rowSums(array(y[,index[j]:(index[j]-1+L[j])],dim=c(n,1)))
#theta<-as.numeric(glm(sumy~x[,],family=poisson,weights = z[,k])$coefficients)
sumy<-rowSums(array(y[sam,index[j]:(index[j]-1+L[j])],dim=c(length(sam),1)))
theta<-as.numeric(glm(sumy~x[sam,],family=poisson,weights = z[sam,k]/previousz[sam,big.cluster])$coefficients)
theta1<-theta
# Main Loop of N-R iterations
metritis<-1
#print(c("dwstou"))
while (sc > nrthreshold&metritis<maxnr ) {
#print(theta)
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
#print(sc)
if (is.nan(sc)==T)
#{theta<-runif(1+tau);sc<-35} else 
{theta<-theta1;sc<-35} else  
if (is.infinite(sc)==T)
#{theta<-runif(1+tau);sc<-35} else
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T){theta<-theta1} #else theta<-runif(1+tau)
#if (sc>gr){theta<-runif(1+tau);print("oops")} else {
#theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#}
}
#print(metritis)
alpha[iter,j,k] <- theta[1]# alpha[j,]'s
beta[iter, j,k,] <- theta[2:(1+tau)]# delta_{j}



}
#############eedw
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
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,k,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}





########################################### liges epanalipseis tou em ###################################################



v<-numeric(K);
yy<-array(data = 0, dim = c(n,q))
for (j in 1:q){
for (l in 1:L[j]){
yy[,j]<-yy[,j] + y[,index[j]-1+l]
}
yy[,j]<-yy[,j]/L[j]
}
st<-numeric(tau + 1)
for (gr in 1:msplit){
############################################################################################################################
#E-Step: Mean allocation vectors   #
############################################################################################################################
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
for(j in 1:q){
#st[1]<-alpha[iter,j,k]
#st[2:(tau+1)]<-beta[iter,j,k,]
#run<-glm(yy[,j]~x,family=poisson, start = st, weights = z[,k])

nrthreshold <- log(10^(-10)); 
#maxnr = 3
sc <- nrthreshold + 1
theta[1] <- alpha[iter,j,k]# alpha[j,]'s
theta[2:(tau+1)] <- beta[iter, j,k,]# beta_{j}
theta1<-theta
# Main Loop of N-R iterations
metritis<-1
#print(c("dwstou"))
while (sc > nrthreshold&metritis<maxnr ) {
#print(theta)
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
#print(sc)
if (is.nan(sc)==T)
#{theta<-runif(1+tau);sc<-35} else 
{theta<-theta1;sc<-35} else  
if (is.infinite(sc)==T)
#{theta<-runif(1+tau);sc<-35} else
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T){theta<-theta1} #else theta<-runif(1+tau)
#if (sc>gr){theta<-runif(1+tau);print("oops")} else {
#theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#}
}
#print(metritis)
alpha[iter,j,k] <- theta[1]# alpha[j,]'s
beta[iter, j,k,] <- theta[2:(1+tau)]# delta_{j}



}
#############eedw
}

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




}


lll<-mylogLikePoisMix(y, condmeans, psim[iter,])



#########################################################################################################################



ll[iter]<-lll[[1]]
print(paste("loglikelihood: ",lll[[1]]))


}

print("end of initialization")

ss<- 1:m
mm<-max(ll)
max.index<- ss[ll==mm]
max.index<-max.index[1]
results<-list(psim[max.index,],alpha[max.index,,], beta[max.index,,,],ll[max.index])
names(results)<-c("psim","alpha,","beta","ll")

}






if (model == 2){

delta <- array(data = 0,dim = c(m,q))

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

if (1<2){

#this is for random selection of observations within the splitted component
#sam<-sample(length(which(clust==big.cluster)),floor(3+runif(1)*(length(which(clust==big.cluster))-3)))
#sam<-index1[sam]
sam<-index1
print(paste("splitting component", big.cluster))

} else{
#this is for kmeans splitting
print(paste("cluster=", big.cluster))
cl<-kmeans(log(rowSums( y[clust==big.cluster,]) +0.5 )/rowSums(array(x[clust==big.cluster,]+0.5,dim=c(length(x[clust==big.cluster,]),tau))),2)
len1<-length(index1[cl$cluster==1])
len2<-length(index1[cl$cluster==2])
print(paste("lengths=", len1,len2))
if (len1*len2==0){
sam<-sample(length(which(clust==big.cluster)),floor(1+length(which(clust==big.cluster))/(1+10*runif(1))))
sam<-index1[sam]
} else{

previousz[index1[cl$cluster==2],big.cluster]<-1
sam<- index1[cl$cluster==1]
}
print(paste("kmeans-splitting component", big.cluster, " n = ", length(sam)))

}

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

nrthreshold <- log(10^(-10));
#maxnr = 3

#######################################################################################################################
#######################################################################################################################
###############################initialization of values################################################################
#######################################################################################################################
#######################################################################################################################
for (j in 1:q){ 
delta[iter,j] <- mean(previous.alpha[j,])
alpha[iter,j,1:(K-1)]<- previous.alpha[j,] - delta[iter,j]
alpha[iter,j,K]<-2*runif(1) - 1
alpha[iter,j,] <- alpha[iter,j,] - mean(alpha[iter,j,]) 
}
beta[iter,,] <- previous.beta

for (k in 1:K) {psim[iter,k]  <- sum(z[,k])/n}

for (j in 1:q){
# initial values
sc <- nrthreshold + 1
theta[1] <- 20*runif(1)# lagrange multiplier 
theta[2:(K+1)] <- alpha[iter,j,]# alpha[j,]'s
theta[K+2] <- delta[iter, j]# delta_{j}
theta[(K+3):(K+2+tau)] <- beta[iter, j,]# beta_{j}
theta1<-theta
#iter1<-0
# Main Loop of N-R iterations
#maxnr=10
metritis = 1
while (sc > nrthreshold&metritis<maxnr ) {
metritis = metritis + 1
grad <- numeric(K + tau + 2)
#grad[1]<-sum(theta[2:(K+1)])
hessian <- array(data = 0, dim = c(K + tau + 2,K + tau + 2))
# partial derivatives wrt Lagrange multiplier and alpha_{j,}
hessian[1,2:(K+1)]<-1;
for (l in 1:L[j]){
for (k in c(K)){
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
{theta<-theta1;sc<-35} else  
if (is.infinite(sc)==T)
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T)theta<-theta1

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
u<-numeric(n);for (t in 1:tau){u<-u + beta[iter,j,t]*x[,t]}
for (l in 1:L[j]){
i<-i+1
ar[,i]<-exp(alpha[iter,j,k] + delta[iter,j] + gamma[j,l] + u)
}}
condmeans[[k]]<-ar
}






###############################################################################################################################

nrthreshold<- -6*log(10)
for (gr in 1:msplit){
############################################################################################################################
#E-Step: Mean allocation vectors   #
############################################################################################################################

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


for (j in 1:q){
# initial values
sc <- nrthreshold + 1
theta[1] <- 20*runif(1)# lagrange multiplier 
theta[2:(K+1)] <- alpha[iter,j,]# alpha[j,]'s
theta[K+2] <- delta[iter, j]# delta_{j}
theta[(K+3):(K+2+tau)] <- beta[iter, j,]# beta_{j}
theta1<-theta
#iter1<-0
# Main Loop of N-R iterations
#maxnr=10
metritis = 1
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
{theta<-theta1;sc<-35} else  
if (is.infinite(sc)==T)
{theta<-theta1;sc<-35} else
theta <- theta - qr.coef(qr(hessian,tol = 1e-300),grad)
#print(theta)
if (is.na(max(theta))==T)theta<-theta1

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




}


lll<-mylogLikePoisMix(y, condmeans, psim[iter,])

##############################################################################################################################






ll[iter]<-lll[[1]]
print(paste("loglikelihood ", lll[[1]]))
}
print("end of initialization")
ss<- 1:m
mm<-max(ll)
max.index<- ss[ll==mm]
max.index<-max.index[1]
print(c("iter selected = ", max.index))
results<-list(psim[max.index,],alpha[max.index,,], beta[max.index,,],delta[max.index,], ll[max.index])
names(results)<-c("psim","alpha,","beta","delta", "ll")}
return(results)}
