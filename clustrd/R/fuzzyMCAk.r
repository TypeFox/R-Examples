fuzzyMCAk <- function(data,nclus=3,ndim=2,nstart=1){
  
  #in case of binary input
  minobs = min(data)
  maxobs = max(data)
  
  data = data.matrix(data)
  n=nrow(data)
  zitem=ncol(data)

  if ((minobs==0) & (maxobs==1)) {
    zncati=apply(data,2,max)+1 
  } else {
    zncati=apply(data,2,max)
  }
  
  zz=disjMake(data)$dZ  
  zncat= sum(zncati)
  oner = matrix(1,n,1)
  muz  = apply(zz,2,mean)
  z= zz-oner %*% muz  
  rr=apply(t(zz) %*% zz,2,sum )
  Dr=diag(rr)
  sqDr=diag(sqrt(rr))
  invsqDr = diag(1/sqrt(rr))
  Pz=z %*% invsqDr
  PPz=t(Pz) %*% Pz
  svdPz = svd(PPz)
  V=svdPz$v
  invsqD=diag((1/sqrt(svdPz$d)))  
  Fm=Pz %*% V %*% invsqD
  F0=Fm[,c(1:ndim)]  
  oldf=1000000
  
  for(b in 1:nstart){
    Fv={}
    Fm= F0
    U0 <- cmeans(F0,nclus,method="cmeans")$membership
    itmax= 100
    ceps= 0.00001
    imp = 100000
    f0 = 1000000
    it = 0 
    
    while((it <= itmax ) && ( imp > ceps ) ){
      it=it+1
      ####################################################################
      ## STEP 1: update F #############################################
      ####################################################################
      U = U0^2
      DU = matrix(0,nrow=n, ncol=n)
      Psi = matrix(0,nrow=n, ncol=n)
      for (k in 1:nclus)
      {
        Ue = U[,k]
        DU = DU + diag(Ue)
        Psi = Psi + (Ue%*%t(Ue)) /  sum(Ue)
      }
      Omega = matrix(0,nrow=n, ncol=n)
      ii = 0
      for (j in 1:zitem)
      {
        i = ii +1
        ii = ii + zncati[j]
        T = zz[,i:ii]
        Omega = Omega + T %*% pseudoinverse(t(T)%*%T)%*%t(T)
      }
      MU = Omega + Psi - DU
      
      muout = irlba(MU,nv=nclus)
      Fm = muout$v[,c(2:(ndim+1))]
      ####################################################################
      ## STEP 2: update Wj and Center ####################################
      ####################################################################
      f1 = 0
      ii = 0
      A = {}
      for (j in 1:zitem)
      {
        i = ii +1
        ii = ii + zncati[j]
        T=zz[,c(i:ii)]
        W = pseudoinverse(t(T)%*%T)%*%t(T)%*%Fm
        A = rbind(A,W)
        f1 = f1 + sum(t(Fm)%*%Fm - t(Fm)%*%T%*%W)
      }
      f2 = 0
      center = {}
      for (k in 1:nclus)
      {
        Ue = U[,k]
        ct = (t(Ue)%*%Fm)/sum(Ue)
        center = rbind(center,ct)
        f2 = f2 + sum(t(Fm)%*%diag(Ue)%*%Fm) - sum(t(Fm)%*%Ue%*%ct)
      }
      
      cl<-cmeans(Fm,nclus,centers=center,method="cmeans")
      U0 = cl$membership
      index = cl$cluster #hard partition
      
      f=f1+f2
      imp=f0-f
      f0=f
      Fv = cbind(Fv,f)
    } # end WHILE
    
    if (f<= oldf){
      oldF=Fm
      oldindex=index
      oldf=f				 
      Uold=U0
      Aold=A
      centerold = center
    }
  } ##end of FOR
  
  PC = sum(t(sum(U0^2)))/n                  # Partition Coefficient
  FPI = 1-(nclus*PC-1)/(nclus-1)                   # Fuzzyness Performance Index
  PE = -sum(sum(U0*log(U0)))/n            # Partition Entropy
  MPE = PE/log(nclus)                          # Modified Partition Entropy
  
  out=list()
  Fm=oldF
  index=oldindex
  f=oldf
  U=Uold
  A=Aold
  center=centerold
  out$obscoord=Fm
  out$attcoord=A
  out$centroid=center
  out$cluID = as.numeric(index) #hard partition
  out$U = U
  out$FPI = FPI
  out$MPE = MPE
  out$criterion=f
  
  out
}  