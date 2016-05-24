covreg.mcmc <-
function(fmean,fcov,data=NULL,R=1,niter=10000,nthin=10,nsave=niter/nthin,verb=T){

## X1 is the covariate for the mean model ##
## X2 is the covariate for the cov model ##

m1=model.frame(fmean,data)
X1=model.matrix(fmean,m1)
m2=model.frame(fcov,data)
X2=model.matrix(fcov,m2)
Y=model.response(m1)

n=dim(Y)[1]
p=dim(Y)[2]
q1=dim(X1)[2]
q2=dim(X2)[2]

#Y=scale(Y,scale=F)
## starting values ##
B2=array(0,dim=c(p,q2,R))
f0=lm(fmean,data)
B1=t(f0$coef)
G=matrix(rnorm(n*R),ncol=R)
A=cov(f0$resid)
iA=solve(A)


####

## priors ##
nu0=p+2
A0=cov(f0$resid)
iA0=solve(A0)
B01=B1
V01=n*solve(t(X1)%*%X1)
iV01=solve(V01)
V02=n*solve(t(X2)%*%X2)
iV02=solve(V02)

####

############################################

b1.save=array(dim=c(p,q1,nsave))
b2.save=array(dim=c(p,q2,R,nsave))
a.save=array(dim=c(p,p,nsave))
##################################################


## main loop ##
for (ns in 1:niter){

## full conditionals ##

## update B2 and G 
for(k in sample(1:R)){
	E<-Y-X1%*%t(B1) 
	#E=Y
    	for(l in (1:R)[-k]) { E<-E-(X2%*% t(B2[,,l]))*G[,l]  } 
    	s2g<-1/( 1+  diag( X2%*%t(B2[,,k])%*%iA%*%B2[,,k]%*%t(X2)) )
    	eg<-s2g*diag(  X2%*%t(B2[,,k])%*%iA%*%t(E) )  
    	G[,k]<-rnorm(n,eg,sqrt(s2g))
     
    	Xg<-diag(G[,k])%*%X2
    	Bn<- (t(E)%*%Xg)%*%solve(( t(Xg)%*%Xg + iV02 ) )
    	B2[,,k]<-rmn(Bn, A , solve( t(Xg)%*%Xg + iV02 )) 
}
  ## 

## update A and B1 conditional on B2 and G
E<-Y
SSB<-matrix(0,p,p) 
for(k in 1:R){ 
	E<- E - (X2%*% t(B2[,,k]))*G[,k] 
    	SSB<-SSB+ B2[,,k]%*%iV02%*%t(B2[,,k])
}
Bn=(t(E)%*%X1+B01%*%iV01)%*%solve(t(X1)%*%X1 + iV01 ) 
B1=rmn(Bn,A,solve(t(X1)%*%X1+iV01))
E=E-X1%*%t(B1)
SSB=SSB+(B1-B01)%*%iV01%*%t(B1-B01)
An<-A0 + t(E)%*%E + SSB 
iA<-rwish( solve(An), nu0 + n + q1+q2*R  ) ; A<-solve(iA)


## output ##

if (ns%%nthin==0){
b1.save[,,ns/nthin]=B1
b2.save[,,,ns/nthin]=B2
a.save[,,ns/nthin]=A
}
if (verb==T & ns%%(niter/100)==0)cat(round(ns/niter*100),"% done",date(), "\n") 
}


return(list(B1.psamp=b1.save,B2.psamp=b2.save,A.psamp=a.save,matrix.mean=X1,matrix.cov=X2))
}
