
PCDHess.wf<-function(DATA,beta,gamma,nu,order,knots){

###########################
#Second partial derivatives
###########################
t=DATA$t
x=DATA$x
z=DATA$z
k=DATA$k
K=DATA$K

p<-ncol(x)
n<-nrow(x)
L<-order+length(knots)-2

IM<-array(NA,dim=c(L,K,n))
for (i in 1:n){
  for (l in 1:L){
    temp0<-Ispline(t[i,1:k[i]],order,knots)
    IM[l,1:k[i],i]<-temp0[l,]
  }
}

EPHI<-function(beta,gamma,nu){
temp1<-c()
for (i in 1:n){
fz<-nu+sum(z[i,1:k[i]])    #fz:numerator
fm<-nu+gamma%*%IM[,k[i],i]*exp(x[i,]%*%beta) #fm:denominator
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
return(ZM)  #fix i and l return lth vector
}

Ephi<-EPHI(beta,gamma,nu);Ez<-EZ(gamma)

A1<-array(NA,dim=c(p,p,n)) #Generate n p by p matrix and add them together
for (i in 1:n){
A1[,,i]<--c((gamma%*%IM[,k[i],i])*exp(x[i,]%*%beta)*Ephi[i])*(x[i,]%*%t(x[i,]))
}
A<-rowSums(A1,dims=2)   #p by p

B<-matrix(NA,p,L)      #p by L
B1<-array(NA,dim=c(n,p,L))
SEz<-c()  
for (l in 1:L){
 for (i in 1:n){
 B1[i,,l]<-IM[l,k[i],i]*exp(x[i,]%*%beta)*Ephi[i]*x[i,]
 }
B[,l]<--colSums(B1[,,l])
SEz[l]<-sum(Ez[l,,],na.rm =TRUE)
}


judg1<-SEz/gamma^2
#if (any(is.nan(judg1))==TRUE) judg1[which(is.nan(judg1)==TRUE)]<-0
if(any(gamma>=0&gamma<=10^(-30))==TRUE) judg1[which(gamma>=0&gamma<=10^(-30))]<-0;
C<-diag(-judg1)
#C<-diag(-SEz/gamma^2)  #L by L diagnoal matrix
D<-n/nu-n*trigamma(nu) #1 by 1

VARA<-rbind(cbind(A,matrix(0,p,1),B),cbind(matrix(0,1,p),matrix(D,1,1),matrix(0,1,L)),cbind(t(B),matrix(0,L,1),C))

#############################################################
#Variance of the first derivative of the augmented likelihood
#############################################################
#variance and covariance needed for computing covariance matrix

Vphi<-c(); Vlog<-c(); COVphi<-c() 
for (i in 1:n){
fz<-nu+sum(z[i,1:k[i]])                      #fz:numerator
fm<-nu+gamma%*%IM[,k[i],i]*exp(x[i,]%*%beta) #fm:denominator
Vphi[i]<-fz/fm^2
Vlog[i]<-trigamma(fz)
COVphi[i]<-fz/fm*(digamma(fz+1)-digamma(fz))
}

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
COVblk<-matrix(0,S,n)    #S by n
for (i in 1:n){
s=1
 for (lu in 1:(L-1)){
  for (ld in (lu+1):L){
  COVz[s,1:k[i],i]<--COVZM[lu,1:k[i],i]*COVZM[ld,1:k[i],i]*z[i,1:k[i]]
  COVblk[s,i]<-IM[lu,k[i],i]*IM[ld,k[i],i]*exp(x[i,]%*%beta)^2*Vphi[i]
  s=s+1
  }
 }
}

E1<-array(NA,dim=c(p,p,n))
G1<-matrix(NA,p,n)
for (i in 1:n){
E1[,,i]<-c(((gamma%*%IM[,k[i],i])*exp(x[i,]%*%beta))^2*Vphi[i])*(x[i,]%*%t(x[i,]))
G1[,i]<-gamma%*%IM[,k[i],i]*exp(x[i,]%*%beta)*(Vphi[i]-COVphi[i])*x[i,]
}
E<-rowSums(E1,dims=2)   #p by p
G<-rowSums(G1)          #p by 1  

F1<-array(NA,dim=c(n,p,L))   # n by p by L
H1<-matrix(NA,L,n)
I1<-matrix(NA,L,n)
SVz<-c()  
F<-matrix(NA,p,L)
for (l in 1:L){
 for (i in 1:n){
 F1[i,,l]<-(gamma%*%IM[,k[i],i])*IM[l,k[i],i]*exp(x[i,]%*%beta)^2*Vphi[i]*x[i,]
 H1[l,i]<-IM[l,k[i],i]*exp(x[i,]%*%beta)*(Vphi[i]-COVphi[i])
 I1[l,i]<-(IM[l,k[i],i]*exp(x[i,]%*%beta))^2*Vphi[i]
 }
F[,l]<-colSums(F1[,,l])
SVz[l]<-sum(VZM[l,,],na.rm =TRUE)
}
H<-as.matrix(rowSums(H1),L,1)

judg2<-SVz/gamma^2
if(any(gamma>=0&gamma<=10^(-30))==TRUE) judg2[which(gamma>=0&gamma<=10^(-30))]<-0;
I<-diag(rowSums(I1)+judg2)

SCOVz<-c()
for (ll in 1:S){
SCOVz[ll]<-sum(COVz[ll,,]) #gamma_l*gamma_l' canceled
}
if (any(is.nan(SCOVz))==TRUE) SCOVz[which(is.nan(SCOVz)==TRUE)]<-0;
J<-rowSums(COVblk)+SCOVz

I[lower.tri(I)]<-J    #combine I and J into one L by L matrix
temp2<-t(I)           #upper triangular
temp3<-I+temp2
diag(temp3)<-diag(temp3)/2 #actual I+J

M<-sum(Vlog)-2*sum(COVphi)+sum(Vphi)

VARB<-rbind(cbind(E,as.matrix(G,p,1),F),cbind(matrix(G,1,p),matrix(M,1,1),t(H)),cbind(t(F),H,temp3))

HESS<-(-VARA-VARB)
return(HESS)
}


