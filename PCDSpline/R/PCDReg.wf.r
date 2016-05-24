
PCDReg.wf<-function(DATA,order,placement,nknot,myknots,binit,ninit,ginit,t.seq,tol){

L=order+nknot-2   #number of basis functions

t=DATA$t
x=DATA$x
z=DATA$z
k=DATA$k
K=DATA$K

#equally-spaced knots
if (placement==TRUE){
t.max<-max(t,na.rm=TRUE)+0.00001
t.min<-min(t,na.rm=TRUE)-0.00001
knots<-seq(t.min,t.max,length.out=nknot)   
}

#specified knots by user
if (placement==FALSE){
knots<-myknots
}

n=nrow(x)
IM<-array(NA,dim=c(L,K,n))
for (i in 1:n){
 for (l in 1:L){
  temp0<-Ispline(t[i,1:k[i]],order,knots)
  IM[l,1:k[i],i]<-temp0[l,]
 }
}


###################
#Three expectations
###################

EPHI<-function(beta,gamma,nu){
temp1<-c()
for (i in 1:n){
fz<-nu+sum(z[i,1:k[i]])   
fm<-nu+gamma%*%IM[,k[i],i]*exp(x[i,]%*%beta)
temp1[i]<-fz/fm
}
return(temp1)
}

EZ<-function(gamma){
ZM<-array(NA,dim=c(L,K,n))
for (i in 1:n){
 for (l in 1:L){
 ZM[l,1,i]<-gamma[l]*IM[l,1,i]*z[i,1]/gamma%*%IM[,1,i]
if (k[i]>1){
 ZM[l,2:k[i],i]<-gamma[l]*(IM[l,2:k[i],i]-IM[l,1:(k[i]-1),i])*z[i,2:k[i]]/(gamma%*%IM[,2:k[i],i]-gamma%*%IM[,1:(k[i]-1),i])
  }
 }
}
return(ZM)  
}

ELOG<-function(beta,gamma,nu){
temp2<-c()
for (i in 1:n){
temp2[i]=digamma(nu+sum(z[i,1:k[i]]))-log(nu+gamma%*%IM[,k[i],i]*exp(x[i,]%*%beta))
}
return(temp2)
}

################
#Three Equations
################

###Beta Function
BETAF<-function(beta,beta0,gamma0,nu0){
Ephi<-EPHI(beta0,gamma0,nu0)
Ez<-EZ(gamma0)
ebeta<-exp(x%*%beta)

blK<-matrix(,n,L)     
SEz<-c()   
for (l in 1:L){
 for (i in 1:n){            
 blK[i,l]<-IM[l,k[i],i]
 }
SEz[l]<-sum(Ez[l,,],na.rm=TRUE)
}
temp3<-(blK%*%diag(SEz)*c(ebeta*Ephi))%*%(diag(1/colSums(blK*c(ebeta*Ephi))))
temp4<-(rowSums(z,na.rm=TRUE)-rowSums(temp3))*x
colSums(temp4)
}

###Gamma function
GAMMAF<-function(beta1,beta0,gamma0,nu0){
ebeta<-exp(x%*%beta1)
Ephi<-EPHI(beta0,gamma0,nu0)
Ez<-EZ(gamma0)

blK<-matrix(,n,L)      
SEz<-c()   
for (l in 1:L){
 for (i in 1:n){             
 blK[i,l]<-IM[l,k[i],i]
 }
SEz[l]<-sum(Ez[l,,],na.rm =TRUE)
}
SEz/colSums(blK*c(ebeta*Ephi))
}

#Nu Function
NUF<-function(beta0,gamma0,nu0,nu){
temp5<-nu*sum(ELOG(beta0,gamma0,nu0))-nu*sum(EPHI(beta0,gamma0,nu0))+n*nu*log(nu)-n*log(gamma(nu))
return(-temp5)
}

###Start Iterations
b0<-binit      #starting value of beta 
g0<-ginit      #starting value of gamma
n0<-ninit
tol<-tol

b1<-nleqslv(b0,BETAF,method="Newton",beta0=b0,gamma0=g0,nu0=n0)$x
g1<-GAMMAF(beta1=b1,beta0=b0,gamma0=g0,nu0=n0)
n1<-optimize(NUF,c(0.1,10),beta0=b0,gamma0=g0,nu0=n0)$minimum

while (max(abs(c(b0,g0,n0)-c(b1,g1,n1)))>tol){
b0<-b1
g0<-g1
n0<-n1
b1<-nleqslv(b0,BETAF,method="Newton",beta0=b0,gamma0=g0,nu0=n0)$x
g1<-GAMMAF(beta1=b1,beta0=b0,gamma0=g0,nu0=n0)
n1<-optimize(NUF,c(0.1,10),beta0=b0,gamma0=g0,nu0=n0)$minimum
}

parest<-c(b1,n1,g1)

vc1<-PCDHess.wf(DATA,b1,g1,n1,order,knots) 
BASIS<-Ispline(x=t.seq,order=order,knots=knots)
bmf<-colSums(g1*BASIS)

flag=is.non.singular.matrix(vc1)

if (any(g1<=10^(-30))==TRUE) {
 pois<-which(g1==0)
 vc2<-vc1[-(length(b1)+length(n1)+pois),-(length(b1)+length(n1)+pois)]
 varcov<-solve(vc2)
 var.bn<-varcov[1:(ncol(x)+1),1:(ncol(x)+1)]
}else{
 varcov<-solve(vc1)  
 var.bn<-varcov[1:(ncol(x)+1),1:(ncol(x)+1)]
}

LL1<-c();LL2<-matrix(,n,K)
for (i in 1:n){
LL1[i]<-log(n1^n1*gamma(n1+sum(z[i,1:k[i]]))/(gamma(n1)*(n1+g1%*%IM[,k[i],i]*exp(x[i,]%*%b1))^(n1+sum(z[i,1:k[i]]))))
#LL2[i,1]<-log(g1%*%IM[,1,i]*exp(x[i,]%*%b1)/factorial(z[i,1]))
LL2[i,1]<-log(g1%*%IM[,1,i]*exp(x[i,]%*%b1))
 if (k[i]>1){
#LL2[i,2:k[i]]<-log((g1%*%IM[,2:k[i],i]-g1%*%IM[,1:(k[i]-1),i])*c(exp(x[i,]%*%b1))/factorial(z[i,2:k[i]]))
LL2[i,2:k[i]]<-log((g1%*%IM[,2:k[i],i]-g1%*%IM[,1:(k[i]-1),i])*c(exp(x[i,]%*%b1)))
 }
}
ll<-sum(LL1)+sum(LL2,na.rm=TRUE)
AIC<-2*length(parest)-2*ll
BIC<-length(parest)*log(n)-2*ll

return(list("beta"=b1,"nu"=n1,"gamma"=g1,"var.bn"=var.bn,"knots"=knots,"bmf"=bmf,"AIC"=AIC,"BIC"=BIC,"flag"=flag))
}
