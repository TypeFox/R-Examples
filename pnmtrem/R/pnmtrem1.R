pnmtrem1 <-
function(covmat1,covmat2,respmat1,respmat2,z,nsubj,nresp,param01,param02,beta0,alpha0,
tol1=0.0001,tol2=0.0001,maxiter1=50,maxiter2=50,tun1=1,tun2=1,x01=0,eps1=10^-10,x02=0,eps2=10^-10,silent=TRUE,delta.print=FALSE,deltastar.print=FALSE){

covmat1<-as.matrix(covmat1)
covmat2<-as.matrix(covmat2)
respmat1<-as.matrix(respmat1)
respmat2<-as.matrix(respmat2)

# Newton-Raphson Algortihm Function (Ilk,O. 2011)
newton<-function(fun,derf,x0,eps){
iter<-0
repeat{
iter<-iter+1
x2<-x0-fun(x0)/derf(x0)
if(abs(x0-x2)<eps) break
x0<-x2
#cat("Iteration:",iter,"Result:",x2,fill=T)
}
x2
}

#20 point Gauss-Hermite Quadratures and Weights
q<-c(.24534,.7374,1.234,1.73853,2.25497,2.7888,3.3478,3.9447,
4.6036,5.3874,-.24534,-.73747,-1.2340,-1.7385,-2.2549,-2.788,
-3.3478,-3.944,-4.603,-5.387)
w<-c(.4622,.286,.109,.0248,.00324,.00023,.0000078,.0000001,
.00000000044,.0000000000002,.4622,.2866,.10901,.0248,.0032,
.00023,.0000078,.00000011,.00000000044,.00000000000023 )
q<-sqrt(2)*q

library(MASS)

#--------------------
# Initial State Model
#--------------------

param0<-param01
param1<-rep(1,length(param0))
bstar0<-rep(0,dim(covmat1)[2])

iter1<-1
while((t(param0-param1)%*%(param0-param1)>tol1) & (iter1<maxiter1)){

param1<-param0

bstar<-param0[1:length(bstar0)]
lstar<-param0[(length(bstar0)+1):(length(bstar0)+(nresp-1))]
lstar<-c(1,lstar)     #fixing lambda1=1
lsig1<-param0[length(param1)]

### derivative of loglikelihood wrt betastar

sum1<-0
for (i in 1:nsubj){
## h_inverse
sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1=sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
hinv<-1/sum2

## the expression right to the h_inverse
sum4<-0
for (m in 1:length(q)){
sum5<-0
sum6<-0
for (n in 1:nresp){
pp1=sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*
(respmat1[(i+(nsubj)*(n-1)),]-pnorm(pp1)))/
(pnorm(pp1)*(1-pnorm(pp1)))*(sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),])
sum5<-sum5+exp2
} ## n
for(t in 1:nresp){
pp1=sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+sum5*exp(sum6)*w[m]
} ## m
sum1<-sum1+hinv*sum4
}## i
derloglikbstar<-sum1   

### derivative of loglikelihood according to lstar

derlogliklstar<-c(rep(0,(nresp-1)))
for (ls in 2:nresp){
sum1<-0
for (i in 1:nsubj){
## h_inverse
sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1=sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
hinv<-1/sum2

## the expression right to the h_inverse
sum4<-0
for (m in 1:length(q)){
n<-ls
pp1=sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*(respmat1[(i+(nsubj)*(n-1)),]-pnorm(pp1)))/
(pnorm(pp1)*(1-pnorm(pp1)))*
((1+(lstar[n]^2)*exp(2*lsig1))^(-1/2)*(lstar[n])*exp(2*lsig1)*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+q[m]*exp(lsig1))
sum6<-0
for(t in 1:nresp){
pp1=sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+exp2*exp(sum6)*w[m]
} ## m
sum1<-sum1+hinv*sum4
}## i
derlogliklstar[ls-1]<-sum1
}##ls



### derivative of loglikelihood wrt logsigma1

sum1<-0
for (i in 1:nsubj){
## h_inverse
sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1<-sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
hinv<-1/sum2

## the expression right to the h_inverse
sum4<-0
for (m in 1:length(q)){
sum5<-0
sum6<-0
for (n in 1:nresp){
pp1<-sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*(respmat1[(i+(nsubj)*(n-1)),]-
pnorm(pp1)))/
(pnorm(pp1)*
(1-pnorm(pp1)))*
((1+(lstar[n]^2)*exp(2*lsig1))^(-1/2)*(lstar[n]^2)*exp(2*lsig1)*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1))
sum5<-sum5+exp2
} ## n
for(t in 1:nresp){
pp1<-sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+sum5*exp(sum6)*w[m]
} ## m

sum1<-sum1+hinv*sum4
}## i
derlogliksig1<-sum1

derlogliktheta<-c(derloglikbstar,derlogliklstar,derlogliksig1)

#maximized loglikelihood
sum13<-0
for (i in 1:nsubj){
## h
sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1<-sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
logh<-log(sum2)
sum13<-sum13+logh
}
maxloglik1<-sum13


## Information Matrix
sum12<-0

for (i in 1:nsubj){

## hinv_betastar

sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1<-sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
hinvbstar<-1/sum2

## der_h_betastar

sum4<-0
for (m in 1:length(q)){
sum5<-0
sum6<-0
for (n in 1:nresp){
pp1<-sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*(respmat1[(i+(nsubj)*(n-1)),]-
pnorm(pp1)))/
(pnorm(pp1)*
(1-pnorm(pp1)))*
(sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),])
sum5<-sum5+exp2
} ## n
for(t in 1:nresp){
pp1<-sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+sum5*exp(sum6)*w[m]
} ## m
dhbstar<-sum4

#der_lstar
dhlstar<-c(rep(0,(nresp-1)))
for (ls in 2:nresp){
sum4<-0
for (m in 1:length(q)){
n<-ls
pp1=sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*(respmat1[(i+(nsubj)*(n-1)),]-pnorm(pp1)))/
(pnorm(pp1)*(1-pnorm(pp1)))*
((1+(lstar[n]^2)*exp(2*lsig1))^(-1/2)*(lstar[n])*exp(2*lsig1)*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+q[m]*exp(lsig1))
sum6<-0
for(t in 1:nresp){
pp1=sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+exp2*exp(sum6)*w[m]
} ## m
dhlstar[ls-1]<-sum4
} ## ls


## der_h_sig1

sum4<-0
for (m in 1:length(q)){
sum5<-0
sum6<-0
for (n in 1:nresp){
pp1<-sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1)
exp2<-(dnorm(pp1)*(respmat1[(i+(nsubj)*(n-1)),]-
pnorm(sqrt(1+(lstar[n]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1))))/
(pnorm(pp1)*
(1-pnorm(pp1)))*
((1+(lstar[n]^2)*exp(2*lsig1))^(-1/2)*(lstar[n]^2)*exp(2*lsig1)*covmat1[(i+(nsubj)*(n-1)),]%*%bstar+lstar[n]*q[m]*exp(lsig1))
sum5<-sum5+exp2
} ## n
for(t in 1:nresp){
pp1<-sqrt(1+(lstar[t]^2)*exp(2*lsig1))*covmat1[(i+(nsubj)*(t-1)),]%*%bstar+lstar[t]*q[m]*exp(lsig1)
exp3<-respmat1[(i+(nsubj)*(t-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(t-1)),])*(log(1-pnorm(pp1)))
sum6<-sum6+exp3
} ## t
sum4<-sum4+sum5*exp(sum6)*w[m]
} ## m
dhsig1<-sum4

hinvtheta<-c(rep(hinvbstar,length(param0)))
dhtheta<-c(dhbstar,dhlstar,dhsig1)
sum12<-sum12+matrix(hinvtheta*dhtheta)%*%t(matrix(hinvtheta*dhtheta))
} ## i
infmat<-sum12

param2<-matrix(param1)+ginv(infmat)%*%matrix(derlogliktheta)/tun1
param0<-as.numeric(param2)

if (silent==FALSE){
print("baseline")
cat("Iteration:",(iter1-1),"\n")
cat("Begins from:",param1,"\n")
cat("Derivative of loglikelihood:",derlogliktheta,"\n")
cat("Maximized loglikelihood:",maxloglik1,"\n")
cat("Goes to:",param0,"\n")
cat("Difference:",t(param0-param1)%*%(param0-param1),"\n")
}

iter1<-iter1+1
}##end of while

if (iter1>=maxiter1){
stop("Algorithm did not converge for the initial state model") 
}else{
rlen<-length(bstar0)
bstar<-param0[1:rlen]
lstar<-c(1,param0[(rlen+1):(rlen+nresp-1)])
lsig1<-param0[length(param0)]
#maximized loglikelihood
sum13<-0
for (i in 1:nsubj){
## h
sum2<-0
for (j in 1:length(q)){
sum3<-0
for (k in 1:nresp){
pp1<-sqrt(1+(lstar[k]^2)*exp(2*lsig1))*(covmat1[(i+(nsubj)*(k-1)),]%*%bstar)+lstar[k]*q[j]*exp(lsig1)
exp1<-respmat1[(i+(nsubj)*(k-1)),]*log(pnorm(pp1))+
(1-respmat1[(i+(nsubj)*(k-1)),])*log(1-pnorm(pp1))
sum3<-sum3+exp1
}## k
sum2<-sum2+exp(sum3)*w[j]
}## j
logh<-log(sum2)
sum13<-sum13+logh
}
maxloglik1<-sum13

result1<-matrix(rep(NA,(rlen+nresp-1+1)*4),ncol=4)
result1[,1]<-t(param0)
result1[,2]<-sqrt(diag(ginv(infmat)))
result1[1:rlen,3]<-result1[1:rlen,1]/result1[1:rlen,2]
result1[(rlen+1):(nrow(result1)-1),3]<-(result1[(rlen+1):(nrow(result1)-1),1]-1)/result1[(rlen+1):(nrow(result1)-1),2]
result1[,4]<-(1-pnorm(abs(result1[,3])))*2
colnames(result1)<-c("Estimate","Std. Error", "Z", "P value")
rownames(result1)<-c(paste("betastar",0:(rlen-1),sep=""),paste("lambdastar",2:nresp,sep=""),"logsigma1")
bsfit<-result1[1:rlen,1] # fitted values of betastars
names(bsfit)<-NULL
}


#-----------------------
# Modeling data for t>=2 
#-----------------------

covmat<-covmat2
mresp<-respmat2
param2<-rep(1,length(param02))
beta<-rep(0,dim(covmat)[2])
alpha0<-matrix(alpha0,ncol=ncol(z),byrow=T)

# number of time and minus 1
nt<-dim(covmat)[1]/(nsubj*nresp)+1
nt2<-nt-1

blen<-length(beta) #length of beta
qlen<-length(q) #length of q and m

# creating delta0 matrices to be filled 
delta01<-matrix(rep(0,nresp*nsubj),ncol=1)           # for t=2
delta02<-matrix(rep(0,((nt2-1)*nresp*nsubj)),ncol=1) # for t>=3

# artifical (all elemets 0) matrix added to the beginning of the
# covariate matrix of t>=2
covmat_add<-matrix(rep(0,nresp*nsubj*blen),ncol=blen)
covmatar<-rbind(covmat_add,covmat)

# obtaining delta01 -- t=2
t<-2
for (m in 1:(nresp*nsubj)){
pp1<-covmat1[m,]%*%bsfit
pp2<-z[(t-1)*nresp*nsubj+m-nresp*nsubj,]
fun<-function(x){
pnorm(covmatar[(t-1)*nresp*nsubj+m,]%*%beta0)-pnorm(x)*(1-pnorm(pp1))-pnorm(x+pp2%*%alpha0[t-1,])*pnorm(pp1)
}#fun
derf<-function(x){
-dnorm(x)*(1-pnorm(pp1))-dnorm(x+pp2%*%alpha0[t-1,])*pnorm(pp1)
}#derf
delta01[m,]<-newton(fun,derf,x01,eps1)
}#m

# obtaining delta02 -- t>=3
for (t in 3:nt){
for (m in 1:(nresp*nsubj)){
pp1<-covmatar[(t-2)*nresp*nsubj+m,]%*%beta0
pp2<-z[(t-1)*nresp*nsubj+m-nresp*nsubj,]
fun<-function(x){
pnorm(covmatar[(t-1)*nresp*nsubj+m,]%*%beta0)-pnorm(x)*(1-pnorm(pp1))-pnorm(x+pp2%*%alpha0[t-1,])*pnorm(pp1)
}#fun
derf<-function(x){
-dnorm(x)*(1-pnorm(pp1))-dnorm(x+pp2%*%alpha0[t-1,])*pnorm(pp1)
}#derf
delta02[(t-3)*nresp*nsubj+m,]<-newton(fun,derf,x01,eps1)
}#m
}#t

# delta0
#delta0<-rbind(delta01,delta02)

# obtaining A & B

# z times alpha0
zalpha0<-matrix(z[,1])
for (t in 1:nt2){
zalpha0[((t-1)*nresp*nsubj):(t*nresp*nsubj),]<-z[((t-1)*nresp*nsubj):(t*nresp*nsubj),]%*%matrix(alpha0[t,])
}

# dividing zalpha0 into t=2 and t>=3
zalpha02<-matrix(zalpha0[1:(nsubj*nresp),],ncol=1)
zalpha03<-matrix(zalpha0[(nsubj*nresp+1):(dim(z)[1]),],ncol=1)

# t=2
cov2<-covmat[1:(nresp*nsubj),]
a1<-cov2*as.numeric(dnorm(cov2%*%beta0))
#a2<-covmat1*as.numeric(pnorm(delta01)*dnorm(covmat1%*%bsfit))
#a3<-covmat1*as.numeric(pnorm(delta01+zalpha02)*dnorm(covmat1%*%bsfit))
ab1<-dnorm(delta01)*(1-pnorm(covmat1%*%bsfit))
ab2<-dnorm(delta01+zalpha02)*pnorm(covmat1%*%bsfit)
b1<--(matrix(z[1:(nsubj*nresp),],ncol=dim(z)[2])*as.numeric(dnorm(delta01+zalpha02)*pnorm(covmat1%*%bsfit)))
#A1<-(a1+a2-a3)/as.numeric(ab1+ab2)
A1<-(a1)/as.numeric(ab1+ab2)
B1<-b1/as.numeric(ab1+ab2)

# t>=3
# covariates at time t-1 (covt1) & covariates at time t (covt)
covt1<-covmat[1:(dim(covmat)[1]-nresp*nsubj),]
covt<-covmat[(nresp*nsubj+1):dim(covmat)[1],]

a1<-covt*as.numeric(dnorm(covt%*%beta0))
a2<-covt1*as.numeric(pnorm(delta02)*dnorm(covt1%*%beta0))
a3<-covt1*as.numeric(pnorm(delta02+zalpha03)*dnorm(covt1%*%beta0))
ab1<-dnorm(delta02)*(1-pnorm(covt1%*%beta0))
ab2<-dnorm(delta02+zalpha03)*pnorm(covt1%*%beta0)
b1<--(matrix(z[(nsubj*nresp+1):(dim(z)[1]),],ncol=dim(z)[2])*as.numeric(dnorm(delta02+zalpha03)*pnorm(covt1%*%beta0)))

A2<-(a1+a2-a3)/as.numeric(ab1+ab2)
B2<-b1/as.numeric(ab1+ab2)

A<-rbind(A1,A2)
B<-rbind(B1,B2)

iter2<-1
while((t(param02-param2)%*%(param02-param2)>tol2) & (iter2<maxiter2)){

beta<-param02[1:blen]
alpha<-matrix(param02[(blen+1):(blen+dim(z)[2]*nt2)],ncol=dim(z)[2],byrow=T)
lambda<-param02[(blen+1+dim(z)[2]*nt2):(blen+dim(z)[2]*nt2+nresp-1)]
lambda<-c(1,lambda)
logsig<-param02[(blen+1+dim(z)[2]*nt2+nresp-1):length(param02)]

param2<-param02

##-----------------------------------
## derivatives of the loglikelihood
##-----------------------------------

### derloglik wrt beta
sum1<-0
for (i in 1:nsubj){
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3

# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))*A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum2<-sum2+hinv*sum5
}#t
sum1<-sum1+sum2
}#i

dloglikbeta<-sum1


### derloglik wrt alpha
dloglikalpha<-matrix(rep(0,nt2*dim(z)[2]),ncol=dim(z)[2])
for (t in 2:nt){
sum1<-0
for (i in 1:nsubj){
sum2<-0
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3


# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))*
(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,])

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum1<-sum1+hinv*sum5

}#i

dloglikalpha[t-1,]<-sum1
}#p

dloglikalpha<-as.vector(t(dloglikalpha))


### derloglik wrt logsigma
dlogliklogsigma<-rep(0,nt2)
for (t in 2:nt){
sum1<-0
for (i in 1:nsubj){
sum2<-0
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3


# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

apart<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
bpart<-B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
zpart<-z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
rpart<-mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
derb<-(1/pp1)*(lambda[k]**2)*exp(2*logsig[t-1])*(pp2+pp3+pp5)+pp6

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum1<-sum1+hinv*sum5

}#i

dlogliklogsigma[t-1]<-sum1
}#p

### derloglik wrt lambda
dlogliklambda<-rep(0,(nresp-1))
for (p in 2:nresp){
sum1<-0
for (i in 1:nsubj){
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3

# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
#sum6<-0
#for (k in 1:nresp){
pp1<-sqrt(1+lambda[p]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[p]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(p-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-(pp1)**(-1)*lambda[p]*exp(2*logsig[t-1])*(pp2+pp3+pp5)+exp(logsig[t-1])*q[n]

sum6<-exp2*derb
#sum6<-sum6+exp2*derb
#}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum2<-sum2+hinv*sum5
}#t
sum1<-sum1+sum2
}#i
dlogliklambda[p-1]<-sum1
}

derloglik2<-c(dloglikbeta,dloglikalpha,dlogliklambda,dlogliklogsigma)

##------------------------
## maximized loglikelihood
##------------------------

sum1<-0
for (i in 1:nsubj){
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
sum2<-sum2+log(sum3)
}#t
sum1<-sum1+sum2
}#i
maxloglik2<-sum1

##------------------
##Information Matrix
##------------------

plen<-length(derloglik2)
infmat<-matrix(rep(0,plen*plen),ncol=plen)
for (i in 1:nsubj){

### beta part
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3

# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))*A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum2<-sum2+hinv*sum5
}#t
infbeta<-sum2

### alpha part
infalpha<-matrix(rep(0,nt2*dim(z)[2]),ncol=dim(z)[2])
for (t in 2:nt){
sum1<-0
#for (i in 1:nsubj){
#sum2<-0

# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3


# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))*
(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,])

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum1<-hinv*sum5

#}#i

infalpha[t-1,]<-sum1
}#p

infalpha<-as.vector(t(infalpha))

### logsigma part
inflogsigma<-rep(0,nt2)
for (t in 2:nt){
sum1<-0
#for (i in 1:nsubj){
#sum2<-0
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3


# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
sum6<-0
for (k in 1:nresp){
pp1<-sqrt(1+lambda[k]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[k]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(k-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

apart<-A[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
bpart<-B[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
zpart<-z[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
rpart<-mresp[(t-2)*nresp*nsubj+(k-1)*nsubj+i,]
derb<-(1/pp1)*(lambda[k]**2)*exp(2*logsig[t-1])*(pp2+pp3+pp5)+pp6

sum6<-sum6+exp2*derb
}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum1<-hinv*sum5

#}#i

inflogsigma[t-1]<-sum1
}#p

### lambda part
inflambda<-rep(0,(nresp-1))
for (p in 2:nresp){
sum1<-0
#for (i in 1:nsubj){
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
hinv<-1/sum3

# expression right to hinv
# derivative part
sum5<-0
for (n in 1:qlen){
#sum6<-0
#for (k in 1:nresp){
pp1<-sqrt(1+lambda[p]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(p-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[p]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(p-1)*nsubj+i,]
exp2<-dnorm(ditj)*(yitj-pnorm(ditj))/(pnorm(ditj)*(1-pnorm(ditj)))

derb<-(pp1)**(-1)*lambda[p]*exp(2*logsig[t-1])*(pp2+pp3+pp5)+exp(logsig[t-1])*q[n]

sum6<-exp2*derb
#sum6<-sum6+exp2*derb
#}#k

# l part
sum7<-0
for (s in 1:nresp){
pp1<-sqrt(1+lambda[s]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(s-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[s]*exp(logsig[t-1])*q[n]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(s-1)*nsubj+i,]
exp4<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum7<-sum7+exp4
}#s
sum5<-sum5+sum6*exp(sum7)*w[n]
}#n
sum2<-sum2+hinv*sum5
}#t
#sum1<-sum1+sum2
#}#i
inflambda[p-1]<-sum2
}

infleft<-matrix(c(infbeta,infalpha,inflambda,inflogsigma),ncol=1)
infright<-t(infleft)
infmat<-infmat+infleft%*%infright
}#i

param02<-matrix(param02,ncol=1)+ginv(infmat)%*%matrix(derloglik2,ncol=1)/tun2
param02<-as.vector(param02)

if (silent==FALSE){
print("t>=2")
cat("Iteration:",(iter2-1),"\n")
cat("Begins from:", param2,"\n")
cat("Derivative of loglikelihood:",derloglik2,"\n")
cat("Maximized loglikelihood:",maxloglik2,"\n")
cat("Goes to:",param02,"\n")
cat("Difference:",t(param02-param2)%*%(param02-param2),"\n")
}

iter2<-iter2+1
}##end of while

if (iter2>=maxiter2){
stop("Algorithm did not converge for t>=2")
}else{

beta<-param02[1:blen]
alpha<-matrix(param02[(blen+1):(blen+dim(z)[2]*nt2)],ncol=dim(z)[2],byrow=T)
lambda<-param02[(blen+1+dim(z)[2]*nt2):(blen+dim(z)[2]*nt2+nresp-1)]
lambda<-c(1,lambda)
logsig<-param02[(blen+1+dim(z)[2]*nt2+nresp-1):length(param02)]

##------------------------
## maximized loglikelihood
##------------------------

sum1<-0
for (i in 1:nsubj){
sum2<-0
for (t in 2:nt){
# hinv
sum3<-0
for (m in 1:qlen){
sum4<-0
for (j in 1:nresp){
pp1<-sqrt(1+lambda[j]**2*exp(2*logsig[t-1]))
pp2<--(A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta0+B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%alpha0[t-1,])
pp3<-A[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]%*%beta
pp4<-(B[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]+
z[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]*as.numeric(mresp[(t-2)*nresp*nsubj+(j-1)*nsubj+i,]))
pp5<-pp4%*%alpha[t-1,]
pp6<-lambda[j]*exp(logsig[t-1])*q[m]

ditj<-pp1*(pp2+pp3+pp5)+pp6

yitj<-mresp[(t-1)*nresp*nsubj+(j-1)*nsubj+i,]
exp1<-yitj*log(pnorm(ditj))+(1-yitj)*log(1-pnorm(ditj))
sum4<-sum4+exp1
}#j
sum3<-sum3+exp(sum4)*w[m]
}#m
sum2<-sum2+log(sum3)
}#t
sum1<-sum1+sum2
}#i
maxloglik2<-sum1

result2<-matrix(rep(NA,(blen+dim(z)[2]*nt2+nresp-1+nt2)*4),ncol=4)
result2[,1]<-t(param02)
result2[,2]<-sqrt(diag(ginv(infmat)))
result2[1:(blen+dim(z)[2]*nt2),3]<-result2[1:(blen+dim(z)[2]*nt2),1]/result2[1:(blen+dim(z)[2]*nt2),2]
result2[((blen+dim(z)[2]*nt2)+1):((blen+dim(z)[2]*nt2)+nresp-1),3]<-(result2[((blen+dim(z)[2]*nt2)+1):((blen+dim(z)[2]*nt2)+nresp-1),1]-1)/
result2[((blen+dim(z)[2]*nt2)+1):((blen+dim(z)[2]*nt2)+nresp-1),2]
result2[,4]<-(1-pnorm(abs(result2[,3])))*2
colnames(result2)<-c("Estimate","Std. Error", "Z", "P value")
rownames(result2)<-c(paste("beta",0:(blen-1),sep=""),paste(rep(paste("alpha",2:nt,sep=""),each=dim(z)[2]),1:ncol(z),sep=""),
paste("lambda",2:nresp,sep=""),paste("logsigma",2:nt,sep=""))
}

# Deltastar_i1j
basefit<-result1[,1]
bsfit<-basefit[1:rlen]
lamsfit<-basefit[(rlen+1):(rlen+nresp-1)]
lamsfit<-c(1,lamsfit)
lsig1fit<-basefit[length(basefit)]

delstar1<-NULL
for (i in 1:nresp){
exp1<-sqrt(1+(lamsfit[i]**2)*exp(2*lsig1fit))
exp2<-covmat1[(nsubj*(i-1)+1):(nsubj*i),]%*%bsfit
delstar1<-rbind(delstar1,exp1*exp2)
}# i


# Deltastar_itj & Delta_itj

tgeq2fit<-result2[,1]
bfit<-tgeq2fit[1:blen]
alpfit<-tgeq2fit[(blen+1):(blen+nt2*ncol(z))]
alpfit<-matrix(alpfit,ncol=ncol(z),byrow=T)
lamfit<-tgeq2fit[(blen+nt2*ncol(z)+1):(blen+nt2*ncol(z)+nresp-1)]
lamfit<-c(1,lamfit)
lsigfit<-tgeq2fit[(blen+nt2*ncol(z)+nresp):(length(tgeq2fit))]

delstart<-NULL
delta<-NULL
for (t in 2:nt){
delstart2<-NULL
delta2<-NULL
A1<-A[((t-2)*nsubj*nresp+1):((t-1)*nsubj*nresp),]
B1<-as.matrix(B[((t-2)*nsubj*nresp+1):((t-1)*nsubj*nresp),])
z1<-as.matrix(z[((t-2)*nsubj*nresp+1):((t-1)*nsubj*nresp),])
lag1resp<-as.matrix(mresp[((t-2)*nsubj*nresp+1):((t-1)*nsubj*nresp),])

for (j in 1:nresp){
exp1<-sqrt(1+(lamfit[j]**2)*exp(2*lsigfit[t-1]))
exp2<-A1[((j-1)*nsubj+1):(j*nsubj),]%*%(bfit-beta0)
exp3<-as.matrix(B1[((j-1)*nsubj+1):(j*nsubj),])%*%(alpfit-alpha0)[(t-1),]
exp4<-(as.matrix(z1[((j-1)*nsubj+1):(j*nsubj),])%*%alpfit[(t-1),])*as.matrix(mresp[((j-1)*nsubj+1):(j*nsubj),])
exp5<-exp1*(exp2+exp3+exp4)
exp6<-exp2+exp3
delstart2<-rbind(delstart2,exp5)
delta2<-rbind(delta2,exp6)
}#j
delstart<-rbind(delstart,delstart2)
delta<-rbind(delta,delta2)
}#t

delta<-delta
delstar<-rbind(delstar1,delstart)


lsig1fit<-exp(lsig1fit)
lsigfit<-exp(lsigfit)

# Empirical Bayes
empbay<-NULL
for (i in 1:nsubj){

fun<-function(x){
sum2<-0
for (t in 1:nt){
delstart<-as.matrix(delstar[((t-1)*nsubj*nresp+1):((t)*nsubj*nresp),])
respt<-as.matrix(mresp[((t-1)*nsubj*nresp+1):((t)*nsubj*nresp),])
sum1<-0
for (j in 1:nresp){
delstaritj<-delstart[((j-1)*nsubj+i),]
respitj<-respt[((j-1)*nsubj+i),]
if (t==1){
ditj<-delstaritj+lamsfit[j]*lsig1fit*x
exp1<-lamsfit[j]*lsig1fit*dnorm(ditj)*(respitj-pnorm(ditj))
exp2<-pnorm(ditj)*(1-pnorm(ditj))
exp3<-exp1/exp2
}else{
ditj<-delstaritj+lamfit[j]*lsigfit[t-1]*x
exp1<-lamfit[j]*lsigfit[t-1]*dnorm(ditj)*(respitj-pnorm(ditj))
exp2<-pnorm(ditj)*(1-pnorm(ditj))
exp3<-exp1/exp2
}
sum1<-sum1+exp3
}#j
sum2<-sum2+sum1
}#t
sum2-x
}#fun

derf<-function(x){
sum2<-0
for (t in 1:nt){
delstart<-as.matrix(delstar[((t-1)*nsubj*nresp+1):((t)*nsubj*nresp),])
respt<-as.matrix(mresp[((t-1)*nsubj*nresp+1):((t)*nsubj*nresp),])
sum1<-0
for (j in 1:nresp){
delstaritj<-delstart[((j-1)*nsubj+i),]
respitj<-respt[((j-1)*nsubj+i),]
if (t==1){
ditj<-delstaritj+lamsfit[j]*lsig1fit*x
exp1<--(lamsfit[j]-lsig1fit)**2
exp2<-ditj*(respitj-pnorm(ditj))+dnorm(ditj)
exp3<-pnorm(ditj)*(1-pnorm(ditj))
exp4<-dnorm(ditj)*(respitj-pnorm(ditj))
exp5<-1-2*pnorm(ditj)
exp6<-exp1*(exp2*exp3+exp4*exp5)/(exp3**2)
}else{
ditj<-delstaritj+lamfit[j]*lsigfit[t-1]*x
exp1<--(lamfit[j]-lsigfit[t-1])**2
exp2<-ditj*(respitj-pnorm(ditj))+dnorm(ditj)
exp3<-pnorm(ditj)*(1-pnorm(ditj))
exp4<-dnorm(ditj)*(respitj-pnorm(ditj))
exp5<-1-2*pnorm(ditj)
exp6<-exp1*(exp2*exp3+exp4*exp5)/(exp3**2)
}
sum1<-sum1+exp6
}#j
sum2<-sum2+sum1
}#t
sum2-1
}#der

empbayi<-newton(fun,derf,x02,eps2)
empbay<-c(empbay,empbayi)
}#i


list1<-list()
list1$title<-"FIRST ORDER PROBIT NORMAL MARGINALIZED TRANSITION RANDOM EFFECTS MODELS"
list1$version<-"Package version 1.0"
list1$date<-date()
list1$title1<-"Baseline Model"
list1$output1<-result1
list1$maxloglik1<-maxloglik1
list1$title2<-"t>=2 Model"
list1$output2<-result2
list1$maxloglik2<-maxloglik2
if (delta.print==TRUE) list1$delta<-c(delta)
if (deltastar.print==TRUE) list1$deltastar<-c(delstar)
list1$empbayes<-as.vector(empbay)
list1

} #pnmtrem1

