
PCDReg.nf<-function(DATA,order,placement,nknot,myknots,binit,ginit,t.seq,tol){

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
if (placement==FALSE){knots<-myknots}
  
n=nrow(x)
  
IM<-array(NA,dim=c(L,K,n))
for (i in 1:n){
  for (l in 1:L){
    temp0<-Ispline(t[i,1:k[i]],order,knots)
    IM[l,1:k[i],i]<-temp0[l,]
  }
}
  
############
#Expectation
############
  
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
  
################
#Two Equations
################

###Beta Function
BETAF<-function(beta,beta0,gamma0){
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
  temp3<-(blK%*%diag(SEz)*c(ebeta))%*%(diag(1/colSums(blK*c(ebeta))))
  temp4<-(rowSums(z,na.rm=TRUE)-rowSums(temp3))*x
  colSums(temp4)
}
  
###Gamma function
GAMMAF<-function(beta1,beta0,gamma0){
  ebeta<-exp(x%*%beta1)
  Ez<-EZ(gamma0)
    
  blK<-matrix(,n,L)       
  SEz<-c()   
  for (l in 1:L){
    for (i in 1:n){              
      blK[i,l]<-IM[l,k[i],i]
    }
    SEz[l]<-sum(Ez[l,,],na.rm =TRUE)
  }
  SEz/colSums(blK*c(ebeta))
}
  
###Start Iterations
b0<-binit      #starting value of beta 
g0<-ginit      #starting value of gamma
tol<-tol
  
b1<-nleqslv(b0,BETAF,method="Newton",beta0=b0,gamma0=g0)$x
g1<-GAMMAF(beta1=b1,beta0=b0,gamma0=g0)
  

while (max(abs(c(b0,g0)-c(b1,g1)))>tol){
  b0<-b1
  g0<-g1
  b1<-nleqslv(b0,BETAF,method="Newton",beta0=b0,gamma0=g0)$x
  g1<-GAMMAF(beta1=b1,beta0=b0,gamma0=g0)
}
  
PCDHess.nf<-function(DATA,beta,gamma,IM,L,E1){
  ###########################
  #Second partial derivatives
  ###########################
x=DATA$x
z=DATA$z
k=DATA$k
K=DATA$K
p<-ncol(x)
n<-nrow(x)
Ez<-E1
    
A1<-array(NA,dim=c(p,p,n)) #Generate n p by p matrix and add them together
for (i in 1:n){
  A1[,,i]<--c((gamma%*%IM[,k[i],i])*exp(x[i,]%*%beta))*(x[i,]%*%t(x[i,]))
}
A<-rowSums(A1,dims=2)   #p by p
    
B<-matrix(NA,p,L)      #p by L
B1<-array(NA,dim=c(n,p,L))
SEz<-c()  
for (l in 1:L){
  for (i in 1:n){
    B1[i,,l]<-IM[l,k[i],i]*exp(x[i,]%*%beta)*x[i,]
    }
  B[,l]<--colSums(B1[,,l])
  SEz[l]<-sum(Ez[l,,],na.rm =TRUE)
}
    
judg1<-SEz/gamma^2
#if (any(is.nan(judg1))==TRUE) judg1[which(is.nan(judg1)==TRUE)]<-0
if(any(gamma>=0&gamma<=10^(-30))==TRUE) judg1[which(gamma>=0&gamma<=10^(-30))]<-0;
C<-diag(-judg1)
    
VARA<-rbind(cbind(A,B),cbind(t(B),C))
    
#############################################################
#Variance of the first derivative of the augmented likelihood
#############################################################
#variance and covariance needed for computing covariance matrix
    
VZM<-array(NA,dim=c(L,K,n))
COVZM<-array(NA,dim=c(L,K,n))
for (i in 1:n){
  for (l in 1:L){
    temp0<-gamma[l]*IM[l,1,i]/gamma%*%IM[,1,i]
    VZM[l,1,i]<-temp0*(1-temp0)*z[i,1]
    COVZM[l,1,i]<-temp0/gamma[l]
    if (k[i]>1){
      temp1<-gamma[l]*(IM[l,2:k[i],i]-IM[l,1:(k[i]-1),i])/(gamma%*%IM[,2:k[i],i]-gamma%*%IM[,1:(k[i]-1),i])
      VZM[l,2:k[i],i]<-temp1*(1-temp1)*z[i,2:k[i]]
      COVZM[l,2:k[i],i]<-temp1/gamma[l]
    }
  }
}
    
S<-choose(L,2)
COVz<-array(0,dim=c(S,K,n))  #S by K by n
for (i in 1:n){
  s=1
  for (lu in 1:(L-1)){
     for (ld in (lu+1):L){
      COVz[s,1:k[i],i]<--COVZM[lu,1:k[i],i]*COVZM[ld,1:k[i],i]*z[i,1:k[i]]
      s=s+1
    }
  }
}
    
D<-matrix(0,p,p)   #p by p: second derivative of Beta
E<-matrix(0,p,L)   #p by L: second partial derivativ of Beta and Gamma
    
SVz<-c()  
for (l in 1:L){
  SVz[l]<-sum(VZM[l,,],na.rm =TRUE)
}
    
judg2<-SVz/gamma^2
if(any(gamma>=0&gamma<=10^(-30))==TRUE) judg2[which(gamma>=0&gamma<=10^(-30))]<-0;
F<-diag(judg2)
    
SCOVz<-c()
for (ll in 1:S){
  SCOVz[ll]<-sum(COVz[ll,,])
}
if (any(is.nan(SCOVz))==TRUE) SCOVz[which(is.nan(SCOVz)==TRUE)]<-0;
G<-SCOVz
    
F[lower.tri(F)]<-G    #combine F and G into one L by L matrix
temp2<-t(F)           #upper triangular
temp3<-F+temp2
diag(temp3)<-diag(temp3)/2 #actual F+G
    
VARB<-rbind(cbind(D,E),cbind(t(E),temp3))
    
return(-VARA-VARB)
}
  
parest<-c(b1,g1)
E1<-EZ(g1)
vc1<-PCDHess.nf(DATA,b1,g1,IM,L,E1)  #before inverse matrix
BASIS<-Ispline(x=t.seq,order=order,knots=knots)
bmf<-colSums(g1*BASIS)
  
flag=is.non.singular.matrix(vc1)
  
if (any(g1<=10^(-30))==TRUE) {
pois<-which(g1<=10^(-30))
vc2<-vc1[-(length(b1)+pois),-(length(b1)+pois)]
var.b<-solve(vc2)[1:ncol(x),1:ncol(x)]
}else{
  var.b<-solve(vc1)[1:ncol(x),1:ncol(x)]
}
P<-ncol(x)

#  if (flag) {var.b=solve(vc1)[1:P,1:P]} else {
#  A=vc1[1:P, 1:P]
#  B=vc1[1:P, (P+1):(P+L)]
#  C=vc1[(P+1):(P+L), 1:P]
#  D=vc1[(P+1):(P+L), (P+1):(P+L)]
#  var.b=ginv(A-B%*%ginv(D)%*%C)
#} 
  
LL1<-c();LL2<-matrix(,n,K)
for (i in 1:n){
LL1[i]<-log(exp(-g1%*%IM[,k[i],i]*exp(x[i,]%*%b1)))
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
return(list("beta"=b1,"gamma"=g1,"var.b"=var.b,"Hess"=vc1,"knots"=knots,"bmf"=bmf,"AIC"=AIC,"BIC"=BIC,"flag"=flag))
}

