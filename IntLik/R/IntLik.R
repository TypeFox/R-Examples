#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.


##MCMC Function

#.First.lib<-function(lib,pkg){
#  library.dynam('ILik',pkg,lib)
#}


ILik<-function(L,prior,start,psiseq,psidim=1,proposal="Normal",iternum=1000){

library(maxLik)

m=length(start)

if(proposal== "Normal"){
prop=0
}else if (proposal== "Gamma"){
prop=1
}else{
stop("Invalid proposal value. Proposal can only be Normal or Gamma. ")
}


##rbern function
rbern=function(num,p){
a=rep(0,num)
for (i in 1:num){
b=runif(1,min=-(1-p),max=p)
if(b>0)
a[i]=1
}
a
}

#Normal Proposal
if(prop==0){
#Proposal density
q=function(lambda, lambdap, lambda_star, Sigma_star){
qf=1
for (i in 1:length(lambda)){
qf=qf*dnorm(lambdap[i], mean=lambda_star[i], sd=sqrt(Sigma_star[i,i]))
}
qf
}

#Jump probability alpha
alphaf=function(lambda,lambdap,psi,lambda_star,Sigma_star){
min(1,L(psi,lambdap)*prior(lambdap,psi)/L(psi,lambda)/prior(lambda,psi)*q(lambdap,lambda,lambda_star,Sigma_star)/q(lambda,lambdap,lambda_star,Sigma_star))
}

#Number of Iteration
J=iternum
Jf=1.1*J
Gf=Jf

#Draw sample
SlambdaG=function(psi,J,lambda_star,Sigma_star){
lambdaGf=matrix(0, nrow=Gf, ncol=length(lambda_star))
lambdai=lambda_star
for (j in 1:Jf){
lambdat=lambdai
for(i in 1:m){
lambdat[i]=rnorm(1,mean=lambda_star[i], sd=sqrt(Sigma_star[i,i]))
}
bern=rbern(1,alphaf(lambdai,lambdat,psi,lambda_star,Sigma_star))
if (bern==1){
lambdai=lambdat
}
lambdaGf[j,]=lambdai
}
lambdaGf=lambdaGf[(0.1/1.1*Gf):Gf,]
lambdaGf
}


#Posterior
Pi_hat=function(lambda_star,Sigma_star,psi){
lambdaJf=matrix(0, nrow=J, ncol=length(lambda_star))
lambdaGf=SlambdaG(psi,J,lambda_star,Sigma_star)
num=0
den=0
a=lambda_star
for (j in 1:J) {
for(i in 1:m){
a[i]=rnorm(1, mean=lambda_star[i], sd = sqrt(Sigma_star[i,i]))
}
lambdaJf[j,]=a
}
for (j in 1:J){
if(!is.null(nrow(lambdaGf))){
lambdax=lambdaGf[j,]
}else{
lambdax=lambdaGf[j]
}
num=num+1/J*alphaf(lambdax,lambda_star,psi,lambda_star,Sigma_star)*q(lambdax,lambda_star,lambda_star,Sigma_star)
den=den+1/J*alphaf(lambda_star,lambdaJf[j,],psi,lambda_star,Sigma_star)
}
num/den
}



#Gamma Proposal
}else if (prop==1){
#Proposal density
q=function(lambda,lambdap,lambda_star,Sigma_star){
qf=1
for (i in 1:length(lambda)){
qf=qf*dgamma(lambdap[i],shape=lambda_star[i]^2/Sigma_star[i,i],scale=Sigma_star[i,i]/lambda_star[i])
}
qf
}

#Proposal density Ratio
qR=function(x1,x2,lambda_star,Sigma_star){
qR=1
for (i in 1:m){
k=lambda_star[i]^2/Sigma_star[i,i]
theta=Sigma_star[i,i]/lambda_star[i]
qR=qR*(x1[i]/x2[i])^(k-1)*exp(-1/theta*(x1-x2))
}
}

#Jump probability alpha
alphaf=function(lambda,lambdap,psi,lambda_star,Sigma_star){
min(1,L(psi,lambdap)*prior(lambda,psi)/L(psi,lambda)/prior(lambda,psi)*qR(lambda,lambdap,lambda_star,Sigma_star))
}


#Draw sample
SlambdaG=function(psi,J,lambda_star,Sigma_star){
Jf=1.1*J
Gf=Jf
lambdaGf=matrix(0, nrow=Gf, ncol=length(lambda_star))
lambdai=lambda_star
for (j in 1:Jf){
lambdat=lambdai
for(i in 1:m){
lambdat[i]=rgamma(1, shape=lambda_star[i]^2/Sigma_star[i,i], scale = Sigma_star[i,i]/lambda_star[i])
}
bern=rbern(1,alphaf(lambdai,lambdat,psi,lambda_star,Sigma_star))
if (bern==1){
lambdai=lambdat
}
lambdaGf[j,]=lambdai
}
lambdaGf=lambdaGf[(0.1/1.1*Gf):Gf,]
lambdaGf
}


#Posterior
Pi_hat=function(lambda_star,Sigma_star,psi){
lambdaJf=matrix(0, nrow=J, ncol=length(lambda_star))
lambdaGf=SlambdaG(psi,J,lambda_star,Sigma_star)
num=0
den=0
a=lambda_star
for (j in 1:J) {
for(i in 1:m){
a[i]= rgamma(1, shape=lambda_star[i]^2/Sigma_star[i,i], scale = Sigma_star[i,i]/lambda_star[i])
}
lambdaJf[j,]=a
}
for (j in 1:J){
if(!is.null(nrow(lambdaGf))){
lambdax=lambdaGf[j,]
}else{
lambdax=lambdaGf[j]
}
num=num+1/J*alphaf(lambdax,lambda_star,psi,lambda_star,Sigma_star)*q(lambdax,lambda_star,lambda_star,Sigma_star)
den=den+1/J*alphaf(lambda_star,lambdaJf[j,],psi,lambda_star,Sigma_star)
}
num/den
}

}







#If psi is a scalar but psiseq is a matrix
if(psidim==1 && !is.null(nrow(psiseq)))
 stop("psidim=1 but psiseq is in matrix format.")

if(!is.null(nrow(psiseq))){
if(ncol(psiseq)!=psidim)
 stop("The psidim and the format of psiseq do not match.")
}


#If psi is a scalar and psiseq is a vector
if(psidim==1){
px=psiseq
py=px
#Timer starts
ptm<- proc.time()
for (i in 1:length(px)){
psix=px[i]
#Log-Likelihood(lambda), given psi
logL2<-function(x){
log(L(psix,x))
}

if(prop==0){
a<-maxLik(logLik=logL2, start=start, method="Newton-Raphson")
}else{
A=diag(length(start))
B=rep(0,length(start))
a<-maxLik(logLik=logL2, start=start,  constraints=list(ineqA=A,ineqB=B))
}


lambda_star=coef(a)
Sigma_star=-1/hessian(a)*diag(length(start))
py[i]=exp(log(L(psix,lambda_star))+log(prior(lambda_star,psix))-log(Pi_hat(lambda_star,Sigma_star,psix)))
}
#Timer ends
ptm=proc.time() - ptm
plot(px,py/max(py),main= "Integrated Likelihood",xlab="psi", ylab="Integrated Likelihood")
py

#if psi is a vector and psiseq is a matrix
}else if(psidim>1 && !is.null(nrow(psiseq))){
px=psiseq
py=rep(0,nrow(px))
#Timer starts
ptm<- proc.time()
for (i in 1:nrow(px)){
psix=px[i,]
#Log-Likelihood(lambda), given psi
logL2<-function(x){
log(L(psix,x))
}


if(prop==0){
a<-maxLik(logLik=logL2, start=start, method="Newton-Raphson")
}else{
A=diag(length(start))
B=rep(0,length(start))
a<-maxLik(logLik=logL2, start=start,  constraints=list(ineqA=A,ineqB=B))
}
lambda_star=coef(a)
Sigma_star=-1/hessian(a)*diag(length(start))
py[i]=exp(log(L(psix,lambda_star))+log(prior(lambda_star,psix))-log(Pi_hat(lambda_star,Sigma_star,psix)))
}
#Timer ends
ptm=proc.time() - ptm
py

#if psi is a vector and psiseq is a vector (one psi)
}else if(psidim>1 && length(psiseq)==psidim){
px=psiseq

#Timer starts
ptm<- proc.time()
psix=px
#Log-Likelihood(lambda), given psi
logL2<-function(x){
log(L(psix,x))
}

if(prop==0){
a<-maxLik(logLik=logL2, start=start, method="Newton-Raphson")
}else{
A=diag(length(start))
B=rep(0,length(start))
a<-maxLik(logLik=logL2, start=start,  constraints=list(ineqA=A,ineqB=B))
}
lambda_star=coef(a)
Sigma_star=-1/hessian(a)*diag(length(start))

py=exp(log(L(psix,lambda_star))+log(prior(lambda_star,psix))-log(Pi_hat(lambda_star,Sigma_star,psix)))

#Timer ends
ptm=proc.time() - ptm

py

}else{
 stop("The psidim and the format of psiseq do not match.")
}

}


