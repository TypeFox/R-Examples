####functions for Optimization on Stiefel 
### 5-22-07

###########################I: evaluate basis function on the measurement points
##
Basis.Evl<-function(M,basis.method,t.v,R.inv=NULL,grid=NULL){
##M: dimension to use; basis.method: which basis to use; one of "sin", "poly", "ns", "spike"
##t.v: points to evaluate on
##R.inv, grid: 
##return: evaluation on the basis functions for the ith observation: M by length(t.v)

 result<-numeric(length(t.v))
 if(basis.method=="sin"){
 result<-apply(matrix(t.v),MARGIN=1,Basis.Sin,M=M)
 }
 if (basis.method=="poly"){
 result<-apply(matrix(t.v),MARGIN=1,Basis.Poly,M=M,R.inv=R.inv)
 }

 if (basis.method=="spike"){
 result<-apply(matrix(t.v),MARGIN=1,Basis.Spike,R.inv=R.inv)
 } 

 if (basis.method=="ns"){
 result<-apply(matrix(t.v),MARGIN=1,Basis.NS,grid=grid,B.orth=R.inv)
 } 

 return(result)
}

##
##Phi and auxiliary quantities for ith observation
Phi.aux<-function(M,basis.method,data.list,i,R.inv=NULL, grid=NULL){
##para: M: dimension to use; basis.method: which basis to use; one of "sin", "poly", "ns"
##data.list: sparse, reguar realizations 
##i: ith observation
##return:
#Phi  ##ith component is: M by m_i
#Psi  ##Psi_i=Phi_i%*%t(Phi_i): M by M
#D    ##D_i=Phi_i%*%Y_i%*%t(Y_i)%*%t(Phi_i): M by M

 data.c<-data.list[[i]]
 t.v<-data.c[[1]][,2]
 y.v<-data.c[[1]][,1]
 phi.c<-Basis.Evl(M,basis.method,t.v,R.inv,grid)
 t.phi.c<-t(phi.c) 
 psi.c<-phi.c%*%t.phi.c
 d.c<-phi.c%*%y.v%*%t(y.v)%*%t.phi.c
 return(list(phi.c,psi.c,d.c))
}

#################################II:calculate F_B, gradien(F)_B, and F_BB(delta), 
################################# as function of B, for the ith observation 
###
##FB.HB: caluclate FB and HB(vectorized) for ith observation

  FB.HB<-function(B,phi.aux,sig,Lambda,index,i){
  ##para: B: M by r orthonormal matrix; phi.aux: axuiliary quantities by Phi.aux
  ##sig: error sd; Lambda: eigenvalues
  ##index: Perm.right.index(M,r) for vec(t(X))
  ##i: ith observation
  ##return: 

  M<-nrow(B)
  r<-ncol(B)
  t.B<-t(B)

  temp<-phi.aux[[i]] 
  psi.c<-temp[[2]]
  d.c<-temp[[3]]
  Q.c<-sig^2*diag(1/Lambda)+t.B%*%psi.c%*%B
  Q.inv<-solve(Q.c)
  
  R<-diag(1,M)-psi.c%*%B%*%Q.inv%*%t.B
  S<-psi.c%*%B%*%Q.inv
  T<-d.c%*%B%*%Q.inv

  U1<-2/sig^2*R%*%psi.c 
  V1<-Q.inv%*%t.B%*%T
  U2<-2/sig^2*R%*%(T%*%t.B%*%psi.c-d.c)
  V2<-Q.inv
  U3<-2*R%*%psi.c
  V3<-Q.inv

  W1<-2/sig^2*S
  X1<-R%*%T
  W2<-2/sig^2*R%*%T
  X2<-S
  W3<-(-2)*S
  X3<-S 

  temp1<-kronecker(t(V1),U1)+kronecker(t(V2),U2)+kronecker(t(V3),U3)
  temp2<-kronecker(t(X1),W1)+kronecker(t(X2),W2)+kronecker(t(X3),W3) 
  temp2<-temp2[,index]
  
  result2<-temp1+temp2  ##vector(H(delta))=result2%*%vector(delta); M*r by M*r

##FB
  result1<-as.vector(-2/sig^2*R%*%T+2*S)   #M*r by 1

##   
  result<-cbind(result1,result2)
  return(result)

 }


###
##sum over individual curve results
FB.HB.all<-function(B,phi.aux,sig,Lambda,n,index){
##n: number of observations
M<-nrow(B)
r<-ncol(B)

temp<-apply(matrix(1:n),MARGIN=1,FB.HB, B=B, phi.aux=phi.aux, sig=sig, Lambda=Lambda,index=index)

result<-apply(temp,MARGIN=1,sum)
result<-matrix(result,M*r, 1+M*r)
return(result)
}

###gradient of F, and FBB(as vectorized in front of vercor (delta)) for all data sum together
GradF.FBB<-function(B,phi.aux,sig,Lambda,n,index){
  M<-nrow(B)
  r<-ncol(B)

  temp<-FB.HB.all(B,phi.aux,sig,Lambda,n,index)

  FB<-matrix(temp[,1],M,r)
  gradF<-FB-B%*%t(FB)%*%B

  H<-temp[,-1]
  temp1<-kronecker(t(B),B)[,index]
  
  FBB<-H-temp1%*%H

  return(list(FB,gradF,FBB))
}

####
###equation for delta
Delta.equ<-function(B,phi.aux,sig,Lambda,n,index,cond.tol=1e+10){
##return: A, C: A%*%vec(delta)=C
 M<-nrow(B)
 r<-ncol(B)  
 t.B<-t(B)

 temp<-GradF.FBB(B,phi.aux,sig,Lambda,n,index)
 FB<-temp[[1]]
 gradF<-temp[[2]]
 FBB<-temp[[3]]
 t.FB<-t(FB)

 I<-(-0.5)*B%*%t.FB
 J<-(-0.5)*t.FB%*%B
 eye.M<-diag(1,M)
 PI<-eye.M-B%*%t.B    

 eye.r<-diag(1,r)
 temp1<-kronecker(eye.r,I)+kronecker(t(J),eye.M)+kronecker(J,PI)
 temp2<-0.5*kronecker(t.FB,B)+0.5*kronecker(t.B,FB)
 temp2<-temp2[,index]

 L1<-FBB+temp1+temp2
 L2<-kronecker(eye.r,t.B)+kronecker(t.B,eye.r)[,index] 

 A<-rbind(L1,L2)
 C<-c(-as.vector(gradF),numeric(r^2)) 

##solve for delta
 svd.A<-svd(A)
 U.A<-svd.A$u
 V.A<-svd.A$v
 D.A<-svd.A$d
 if(min(D.A)>1e-3){
 cond.A<-max(D.A)/min(D.A)
 }else{
 cond.A<-2*cond.tol
 }

##check the condition number of A, if it is too big, 
##then stop and return conditional number
 if(cond.A>cond.tol){
  return(list(cond.A))
 } 
 delta<-V.A%*%diag(1/D.A)%*%t(U.A)%*%C
 delta<-matrix(delta,M,r)
 
 return(list(cond.A,A,C,delta,gradF))
} 

####
###update along geodesics with direction delta, with step length sl (0 <=sl <=1)
Update.geo<-function(B, delta,sl=1){
##para: B--current point on Stiefel, M by r, orthornormal
##delta: current direction; tangent vector
##sl: step length
##return: the updated B
 M<-nrow(B)
 r<-ncol(B)  
 t.B<-t(B)
 result<-matrix(0,M,r)
 
 temp<-(diag(1,M)-B%*%t.B)%*%delta
 temp.QR<-qr(temp)
 Q<-qr.Q(temp.QR)
 R<-qr.R(temp.QR)

 A<-t.B%*%delta
 X1<-rbind(A,R)
 X2<-rbind(-t(R),matrix(0,r,r))
 X<-cbind(X1,X2) 
 
 exp.x<-Expo.skew(sl,X)
 
 M.x<-exp.x[1:r,1:r]
 N.x<-exp.x[(r+1):(2*r),1:r]

 result<-B%*%M.x+Q%*%N.x
 return(result) 
}


###############################III: update Lambda and sig^2
###
##calculate gradient and Hessian for log(Lambda) and log(sig^2) for the ith observation
F.Lambda<-function(B,phi.aux,sig,Lambda,data.list,i){
  
  r<-ncol(B)
  t.B<-t(B)
  
  temp<-phi.aux[[i]]
  phi<-temp[[1]]
  t.phi<-t(phi)
  m<-ncol(phi)
  psi.c<-temp[[2]]
  d.c<-temp[[3]]
  Q.c<-sig^2*diag(1/Lambda)+t.B%*%psi.c%*%B
  Q.inv<-solve(Q.c)
  
  P.inv<-(diag(1,m)-t.phi%*%B%*%Q.inv%*%t.B%*%phi)/sig^2
  P.inv2<-P.inv%*%P.inv
  P.inv3<-P.inv2%*%P.inv
  
  data.temp<-data.list[[i]]
  y<-data.temp[[1]][,1]
  t.y<-t(y)

  grad<-numeric(r+1) 
  Hess<-matrix(0,r+1,r+1) 
 
  temp1.0<-(-sig^2)*(t.y%*%P.inv2%*%y)
  temp2.0<-sig^2*sum(diag(P.inv))
  grad[1]<-temp1.0+temp2.0

  temp1.00<-temp1.0+2*sig^4*(t.y%*%P.inv3%*%y)
  temp2.00<-temp2.0-sig^4*sum(diag(P.inv2))
  Hess[1,1]<-temp1.00+temp2.00

  for(k in 1:r){
   bk<-B[,k]
   t.bk<-t(bk)
   temp1.k<-(-Lambda[k])*(t.y%*%P.inv%*%t.phi%*%bk)^2
   temp2.k<-Lambda[k]*t.bk%*%phi%*%P.inv%*%t.phi%*%bk    
   grad[k+1]<-temp1.k+temp2.k
    
   temp1.0k<-2*sig^2*Lambda[k]*(t.y%*%P.inv%*%t.phi%*%bk)*(t.y%*%P.inv2%*%t.phi%*%bk)   
   lk<-t.bk%*%phi%*%P.inv 
   temp2.0k<-(-2*sig^2*Lambda[k])*lk%*%t(lk)
   Hess[1,k+1]<-temp1.0k+temp2.0k
   Hess[k+1,1]<-Hess[1,k+1]
   
   temp1.kk<-temp1.k+2*Lambda[k]^2*(t.y%*%P.inv%*%t.phi%*%bk)^2*(t.bk%*%phi%*%P.inv%*%t.phi%*%bk)
   temp2.kk<-temp2.k-Lambda[k]^2*(t.bk%*%phi%*%P.inv%*%t.phi%*%bk)^2
   Hess[k+1,k+1]<-temp1.kk+temp2.kk 

   for(l in 1:k){
    if(l!=k){      
     bl<-B[,l] 
     temp1.kl<-2*Lambda[k]*Lambda[l]*(t.y%*%P.inv%*%t.phi%*%bk)*(t.y%*%P.inv%*%t.phi%*%bl)*(t(bl)%*%phi%*%P.inv%*%t.phi%*%bk)  
     temp2.kl<-(-Lambda[k]*Lambda[l])*(t.bk%*%phi%*%P.inv%*%t.phi%*%bl)^2
     Hess[k+1,l+1]<-temp1.kl+temp2.kl
     Hess[l+1,k+1]<-Hess[k+1,l+1]
    }   
   }
  }
  
 return(cbind(grad,Hess))
}

####
F.Lambda.all<-function(B,phi.aux,sig,Lambda,data.list,n){
 r<-ncol(B)   
 temp<-apply(matrix(1:n),MARGIN=1,F.Lambda,B=B,phi.aux=phi.aux,sig=sig,Lambda=Lambda, data.list=data.list)
 temp<-apply(temp,1,sum)
 temp<-matrix(temp,r+1,r+2) 
 gradF<-temp[,1]
 Hess<-temp[,-1] 
 return(list(gradF,Hess))
}

####
##update  Lambda and sig
Update.Lambda<-function(B,phi.aux,sig,Lambda,data.list,n,epsi=1e-6){

 ksi0<-log(sig^2) 
 ksi<-log(Lambda)
 
 temp<-F.Lambda.all(B,phi.aux,sig,Lambda,data.list,n)
 gradF<-temp[[1]]
 Hess<-temp[[2]]

  if(sig^2<1e-6){##not update sig 
   gradF.c<-gradF[-1]
   Hess.c<-Hess[-1,-1]
  
   eigen.H<-eigen(Hess.c,symmetric=TRUE)
   U.H<-eigen.H$vectors
   D.H<-eigen.H$values
   H.inv<-U.H%*%diag(1/D.H)%*%t(U.H)
  
   ksi.c<-ksi-H.inv%*%gradF.c 
   Lambda.c<-exp(ksi.c)
   sig.c<-sig

  }else{   
   eigen.H<-eigen(Hess,symmetric=TRUE)
   U.H<-eigen.H$vectors
   D.H<-eigen.H$values
   H.inv<-U.H%*%diag(1/D.H)%*%t(U.H)
  
   temp<-c(ksi0,ksi)-H.inv%*%gradF
   ksi.c<-temp[-1]
   ksi0.c<-temp[1]
   Lambda.c<-exp(ksi.c)
   sig.c<-exp(ksi0.c/2)
  }

 return(list(sig.c,Lambda.c,gradF))
}



###############################IV: Newton's method: iteration and updates
###newton's method with stop rule: either step>max.step, or gradB, gradL<tol;
####or condition> cond.tol
Newton.New<-function(B.c, phi.aux, sig.c, Lambda.c,data.list,n,sl.v,max.step=50,tol=1e-3,cond.tol=1e+10){
 M<-nrow(B.c)
 r<-ncol(B.c)
 index.t<-Perm.right.index(M,r)

##results for return
 like<-numeric(max.step)   ##is this decreasing? not always
  B.result<-array(0,dim=c(M,r,max.step))
  lam.result<-matrix(0,max.step,r)
  sig.result<-numeric(max.step)
  gradB.result<-array(0,dim=c(M,r,max.step))  ##should be zero when converge;
  gradL.result<-matrix(0,max.step,r+1)

##initials 
 error.c<-0
 i<-0
 tol.c<-1e+10
 cond<-(-1)

##loops
 while(i<max.step && tol.c>=tol){
 i<-i+1
 #print(paste("begin iter",i))

 B.result[,,i]<-B.c 
 lam.result[i,]<-Lambda.c 
 sig.result[i]<-sig.c

##calculate likelihood 
 like.c<-try(loglike.all(B.c,phi.aux,sig.c,Lambda.c,data.list,n))
 error.c<-inherits(like.c, "try-error")
  if (error.c){
  like[i]<-(-99)
   break()
  }
  like[i]<-like.c 

 if(i==1){
 gradL.c<-try(F.Lambda.all(B.c,phi.aux,sig.c,Lambda.c,data.list,n)[[1]])
 error.c<-inherits(gradL.c, "try-error")
  if (error.c){
   gradL.result[i,]<-(-99)
   break()
   }
  gradL.result[i,]<-gradL.c
  }


##update B: gradB.result[,,i]
 temp.B<-try(Delta.equ(B.c,phi.aux,sig.c,Lambda.c,n,index.t,cond.tol))
 error.c<-inherits(temp.B, "try-error")
  if (error.c)
   break()

 cond<-temp.B[[1]]
 if(cond>cond.tol)
   break()
 
 gradB.result[,,i]<-temp.B[[5]]
 delta.c<-temp.B[[4]]

 sl.c<-sl.v[i]
 B.up<- try(Update.geo(B.c,delta.c,sl.c))
 error.c<-inherits(B.up, "try-error")
  if (error.c)
   break()

  B.c<-B.up

##update lambda and sig
 temp.Lam<-try(Update.Lambda(B.c,phi.aux,sig.c,Lambda.c,data.list,n,epsi=1e-6))
  error.c<-inherits(temp.Lam, "try-error")
  if (error.c)
   break()

 sig.c<-temp.Lam[[1]]
 Lambda.c<-temp.Lam[[2]] 
 
 if(i<max.step){  
 gradL.result[i+1,]<-temp.Lam[[3]]
 }

## check the abs of gradients 
 tol.c<-max(abs(gradL.result[i,]),abs(gradB.result[,,i]))

##print to screen
 #print(c(paste("end iter",i),round(like[i],3)))
 }

##return
 result<-list(like,B.result,lam.result,sig.result,gradB.result,gradL.result,i,cond,error.c)
 return(result)
}


##################### V: -2loglikehood and CV 
## -2loglikehood  for ith observation 
loglike<-function(B,phi.aux,sig,Lambda,data.list,i){
  M<-nrow(B)
  r<-ncol(B)
  t.B<-t(B)

  temp<-phi.aux[[i]]
  phi<-temp[[1]]
  psi<-temp[[2]]
  Q<-sig^2*diag(1/Lambda)+t.B%*%psi%*%B
  Q.inv<-solve(Q)
  m<-ncol(phi)

  data.temp<-data.list[[i]]
  y<-data.temp[[1]][,1]

  P.inv<-(diag(1,m)-t(phi)%*%B%*%Q.inv%*%t.B%*%phi)/sig^2
  log.det.P<-2*(m-r)*log(sig)+sum(log(Lambda))+log(det(Q))
     
  result<-sum(diag(P.inv%*%y%*%t(y)))+log.det.P+m*log(2*3.1415926)
  return(result) 
}

### the mean -2*loglikelihood over all observations
loglike.all<-function(B,phi.aux,sig,Lambda,data.list,n){
 temp<-apply(matrix(1:n), MARGIN=1, loglike, B=B,phi.aux=phi.aux, sig=sig, Lambda=Lambda, data.list=data.list) 
 result<-mean(temp)
 return(result)
}


###get -2*loglikelihood for ith obs when covariance surface is given on a fine grid
### 
loglike.cov<-function(covmatrix,sig,data.list,i){
  data.temp<-data.list[[i]]
  y<-data.temp[[1]][,1]
  t<-data.temp[[1]][,2]
  m<-length(y)
  timeindex<-floor(t*ncol(covmatrix))+1
 
  P<-sig^2*diag(1,m)+covmatrix[timeindex,timeindex]
  
  temp1<-solve(P,y)%*%t(y)
  temp1<-sum(diag(temp1))
  temp2<-log(det(P))

  result<-temp1+temp2+m*log(2*3.1415926)
  return(result)
}

### mean(-2*loglikelihood)
loglike.cov.all<-function(covmatrix,sig,data.list,n){
 temp<-apply(matrix(1:n), MARGIN=1,loglike.cov,covmatrix=covmatrix,sig=sig,data.list=data.list) 
 result<-mean(temp)
 return(result)
}


####################### cross validation scores
##
CV.B<-function(B,phi.aux,sig,Lambda,n){
###para: B, sig, Lambda: current parameters
##phi.aux: data; n--sample size
##return CV score: the stiefel part; averaged over all obs

M<-nrow(B)
r<-ncol(B)
t.B<-t(B)
index.t<-Perm.right.index(M,r)

##(i) H^{-1}: 
 temp<-apply(matrix(1:n),MARGIN=1,FB.HB,B=B,phi.aux=phi.aux,sig=sig,Lambda=Lambda,index=index.t)
 
 temp.all<-apply(temp,MARGIN=1,sum)
 temp.all<-matrix(temp.all,M*r, 1+M*r)
 FB.all<-matrix(temp.all[,1],M,r)
 gradF.all<-FB.all-B%*%t(FB.all)%*%B
 H.all<-temp.all[,-1]
 temp1<-kronecker(t.B,B)[,index.t]
 FBB.all<-H.all-temp1%*%H.all
 
 I<-(-0.5)*B%*%t(FB.all)
 J<-(-0.5)*t(FB.all)%*%B
 PI<-diag(1,M)-B%*%t.B  
 eye.M<-diag(1,M)
 PI.half<-eye.M-0.5*B%*%t.B   
  
 eye.r<-diag(1,r)
 temp1<-kronecker(eye.r,I)+kronecker(t(J),eye.M)+kronecker(J,PI)
 temp2<-0.5*kronecker(t(FB.all),B)+0.5*kronecker(t.B,FB.all)
 temp2<-temp2[,index.t]

 L1<-FBB.all+temp1+temp2
 L2<-kronecker(eye.r,t.B)+kronecker(t.B,eye.r)[,index.t]

 A<-rbind(L1,L2)
 svd.A<-svd(A)
 U.A<-svd.A$u
 V.A<-svd.A$v
 D.A<-svd.A$d
 temp.A<-V.A%*%diag(1/D.A)%*%t(U.A)   #H^{-1}

##(ii):
partII<-numeric(n)               ##g_c(gradF_i,H^{-1}(gradF_i))
partIII<-numeric(n)              ##H_i(H^{-1}(gradF_i),H^{-1}(gradF_i))

 for (i in 1:n){
 cur<-temp[,i]
 FB.c<-matrix(cur[1:(M*r)],M,r)
 HB.c<-matrix(cur[-(1:(M*r))],M*r,M*r)
 gradF.c<-FB.c-B%*%t(FB.c)%*%B           
  
 ##delta_i=H^{-1}(gradF_i)
 C.c<-c(as.vector(gradF.c),numeric(r^2))  
 delta.c<-temp.A%*%C.c
 delta.c<-matrix(delta.c,M,r)  

 ##FBB_i(H^{-1}(gradF_i),H^{-1}(gradF_i))
 FBB.temp<-matrix(HB.c%*%as.vector(delta.c),M,r)%*%t(delta.c)
 FBB.c<-sum(diag(FBB.temp))  
 
 ##partIII: H_i(delta_i,delta_i)
 temp1<-(t(FB.c)%*%delta.c%*%t.B+t.B%*%delta.c%*%t(FB.c))%*%delta.c
 temp2<-(t.B%*%FB.c+t(FB.c)%*%B)%*%t(delta.c)%*%PI%*%delta.c
 partIII[i]<-1.5*(FBB.c+0.5*sum(diag(temp1))-0.5*sum(diag(temp2)))   

##partII: g_c(gradF_i,delta.c)
 temp.2<-t(gradF.c)%*%PI.half%*%delta.c
 partII[i]<-sum(diag(temp.2))
 }

result<-mean(partII+partIII)
return(result)
}

####
CV.Lambda<-function(B,phi.aux,sig,Lambda,data.list,n){
###para: B, sig, Lambda: current parameters
##phi.aux: data; n--sample size
##return:the CV part on for c(log(Lambda),log(sig^2)); averaged over all obs

 r<-ncol(B)   
 temp<-apply(matrix(1:n),MARGIN=1,F.Lambda,B=B,phi.aux=phi.aux,sig=sig,Lambda=Lambda, data.list=data.list)
 temp.all<-apply(temp,1,sum)
 temp.all<-matrix(temp.all,r+1,r+2) 
 Hess.all<-temp.all[,-1] 

  eigen.H<-eigen(Hess.all,symmetric=TRUE)
  U.H<-eigen.H$vectors
  D.H<-eigen.H$values
  Hess.inv<-U.H%*%diag(1/D.H)%*%t(U.H) 
 
 partII<-numeric(n)
 partIII<-numeric(n)

  for(i in 1:n){
  cur<-temp[,i]
  cur<-matrix(cur,r+1,r+2)
  gradF.c<-cur[,1]
  Hess.c<-cur[,-1]
  
##delta_i=H^{-1}(gradF_i)
  delta.c<-Hess.inv%*%gradF.c

##partII: <gradF_i, delta_i>
  partII[i]<-gradF.c%*%delta.c

##partIII: t(delta.c)%*%Hess.c%*%delta.c
  partIII[i]<-1.5*(t(delta.c)%*%Hess.c%*%delta.c)  
 }

 result<-mean(partII+partIII) 
 return(result)
}


###cross-validation score: averaged over all obs
CV<-function(B,phi.aux,sig,Lambda,data.list,n){

##partI: -2loglikelihood
part.like<-loglike.all(B,phi.aux,sig,Lambda,data.list,n)

##Part for B
part.B<-CV.B(B,phi.aux,sig,Lambda,n)

##part for (log(Lambda),log(sig^2))
part.Lam<-CV.Lambda(B,phi.aux,sig,Lambda,data.list,n)

##
result<-part.like+part.B+part.Lam
return(result)
}

#################################### VII: initial values
Initial<-function(r,ini.method,data.list,n,nmax,grid.l,grids,M.EM=25,iter.num=50,basis.EM="ns",sig.EM){ 

 band<-NULL
 if (ini.method=="loc"){
 temp.loc<-LocLin.Ini(data.list,n,nmax,grid.l,grids,r)
 sig2hat<-temp.loc[[1]]
 eigenf<-temp.loc[[2]]
 eigenv<-temp.loc[[3]]
 band<-temp.loc[[4]] 
 }  

 if(ini.method=="EM"){
 temp.EM<-EM(data.list,n,nmax,grids,M.EM,iter.num,r,basis.EM,sig.EM)
 eigenf<-temp.EM[[1]]*sqrt(length(grids))
 eigenv<-temp.EM[[2]]
 sig2hat<-temp.EM[[3]]^2 
 }
  
 if(ini.method=="Ken"){   ##need to be defined
 }

 covmatrix<-eigenf%*%diag(eigenv)%*%t(eigenf)
 if(sig2hat>0){
  like<-loglike.cov.all(covmatrix,sqrt(sig2hat),data.list,n)
 }else{
 like<-(-99)
 }

result<-list(sig2hat,covmatrix,eigenf,eigenv,like,band)
return(result)
}

###################################auxiliary functions 
###
##permutation of an m*n by m*n matrix, by multiple P_m,n on the right
Perm.right<-function(X,m,n){
##para: X: m*n by m*n matrix
##return: column permutation of X by the rule of P_m,n

index<-matrix(1:(m*n),n,m)
index<-as.vector(t(index))
result<-X[,index]
return(result)
}

Perm.right.index<-function(m,n){
##para: X: m*n by m*n matrix
##return: column permutation of X by the rule of P_m,n

index<-matrix(1:(m*n),n,m)
index<-as.vector(t(index))
return(index)
}

###
##exp(tX) when X is skew-symmetric
Expo.skew<-function(t,X){
 svd.x<-svd(X)
 U.x<-svd.x$u
 D.x<-svd.x$d
 V.x<-svd.x$v

 cos.x<-cos(t*D.x)
 sin.x<-sin(t*D.x)

 result<-U.x%*%diag(cos.x)%*%t(U.x)+U.x%*%diag(sin.x)%*%t(V.x)
 return(result)
}

###
###detect "outlier" singular values of H, avoid singularity, 
###improve stability in computation of newton's method as well as cross validation score
OutDet<-function(sin.H,thre){
##para: sin.H: singular value of H in log scale.
diff.H<-diff(sin.H)
med.H<-median(diff.H)
temp<-abs(diff.H-med.H)
sd.H<-median(temp)
index<-round(length(sin.H)/2):(length(sin.H)-1)
result<-index[temp[index]>thre*sd.H]
if(length(result)==0){
result<-length(sin.H)
}else{
result<-result[1]
}
return(result)  ##return the smallest index falling out of the range 
}




