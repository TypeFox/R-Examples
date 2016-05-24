MCAk <- function(data,nclus=3,ndim=2,nstart=100,smartStart=F,seed=1234){
  data = data.matrix(data)
  
  #in case of binary input
  minobs = min(data)
  maxobs = max(data)
  
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
  M=z
  Pz=z %*% invsqDr
  PPz=t(Pz) %*% Pz
  svdPz = svd(PPz)
  V=svdPz$u
  V=svdPz$v
  D=diag(svdPz$d)
  invsqD=diag((1/sqrt(svdPz$d)))  
  Fm=Pz %*% V %*% invsqD
  F0=Fm[,1:ndim]
  
  oldf=1000000
  
  for(b in 1:nstart){
    Fv={}
    Fm= F0
    
    set.seed(seed+b)
    
    index0=matrix(ceiling(runif(n)*nclus),n,1)
    
    if(smartStart==T){
      index0=kmeans(Fm,nclus,nstart=100)$cluster
    }
    
    U=matrix(0,n,nclus)
    for(j in 1:n){
      U[j,index0[j]]=1
    } 
    
    center= pseudoinverse(t(U) %*% U) %*% t(U) %*% Fm 
    itmax= 100
    it = 0
    ceps= 0.0001
    imp = 100000
    f0 = 100000
    
    itwhile=0 
    
    while((it <= itmax ) && ( imp > ceps ) ){
      it=it+1
      itwhile=itwhile+1
      
      ####################################################################
      ## STEP 1: update of U #############################################
      ####################################################################
      outK = kmeans(Fm,centers=center)
      center=outK$centers
      index = outK$cluster
      U=matrix(0,n,nclus)
      for(i in 1:n){
        U[i,index[i]]=1
      }
      U0=U-oner %*% apply(U,2,mean)
      uu=apply(t(U)%*% U ,2,sum)
      Dru=diag(c(rr,uu))
      sqDru=diag(sqrt(c(rr,uu)))
      invsqDru=diag(1/sqrt(c(rr,uu)))
      
      ####################################################################
      ## STEP 2: update of Fm and Wj ######################################
      ####################################################################
      
      MU=cbind(M,U0)
      Pzu=MU %*% invsqDru
      wzero=(which(Pzu=="NaN",arr.ind=T))
      Pzu[wzero]=0
      PPzu=t(Pzu) %*% Pzu
      
      svdPzu = svd(PPzu)
      Q=svdPzu$u
      
      D=diag(svdPzu$d)
      sqD=diag(sqrt(svdPzu$d))
      invsqD = diag(1/sqrt(svdPzu$d)) 
      
      Fm= Pzu %*% Q %*% invsqD
      Fm= Fm[,1:ndim]
      ft1=0
      k=1
      kk=0
      
      kk=kk+zncati[1]
      Tm=z[,k:kk]
      
      W=pseudoinverse(t(Tm)%*% Tm)%*% t(Tm) %*% Fm
      A=W
      ft1=ft1+sum(diag(t(Fm) %*% Fm-t(Fm) %*% Tm %*% W))
      k=kk+1
      
      for(j in 2:zitem){
        kk=kk+zncati[j]
        Tm=z[,k:kk]
        W=pseudoinverse(t(Tm)%*% Tm)%*% t(Tm) %*% Fm
        A=rbind(A,W)
        ft1=ft1 + sum(diag(t(Fm) %*% Fm-t(Fm) %*% Tm %*% W)) ## MCA criterion
        k=kk+1
      }
      
      ft2=sum(diag(t(Fm) %*% Fm)-t(Fm) %*% U %*% center) ## cluster criterion
      f=ft1+ft2
      imp=f0-f
      f0=f
      Fv = cbind(Fv,f)
    } # end WHILE
    
    if (f<= oldf){
      oldF=Fm
      oldindex=index
      oldf=f				
      Uold=U				
      Aold=A
      centerold = center
    }
  } ##end of FOR
  
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
  out$cluID=as.numeric(index)
  out$criterion=f
  out
}  