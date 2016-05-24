groupals <- function(data,nclus=3,ndim=2,nstart=100,smartStart=F,seed=1234){
  data = data.matrix(data)

  #in case of binary input
  minobs = min(data)
  maxobs = max(data)
  
  outDisj=disjMake(data)
  if ((minobs==0) & (maxobs==1)) {
    zncati=apply(data,2,max)+1 
  } else {
    zncati=apply(data,2,max)
  }
  zncat=sum(zncati)

  q=ncol(data)
  n=nrow(data)
  zz=outDisj$dZ
  
  onen=matrix(1,nrow=n,ncol=1)
  mu_zz=apply(zz,2,mean)
  z=zz - onen %*% mu_zz
  rr=diag(t(zz) %*% zz)
  Dr=diag(rr)
  invDr=diag(1/sqrt(rr))
  Pz=z %*% invDr
  svdRes=svd(t(Pz)%*% Pz)
  V=svdRes$v
  D=diag(svdRes$d)
  invSqD=diag(1/sqrt(svdRes$d))
  U=svdRes$u
  Fv=Pz%*%V %*% invSqD
  F0=Fv[,1:ndim]
  
  oldf = 1000000                   
  
  for (b in 1:nstart){
    Fv=F0
    set.seed(seed+b)
    index0=matrix(ceiling(runif(n)*nclus),n,1)
    
    if(smartStart==T){
      index0=kmeans(Fv,nclus,nstart=100)$cluster
    }
    
    U=matrix(0,n,nclus)
    for(i in 1:n){
      U[i,index0[i]]=1
    }
    
    A=matrix(0,nrow=zncat, ncol=ndim)
    itmax=100
    it=0
    ceps = 0.00001           
    imp = 100000            
    f0 = 1000000
    
    while ((it <= itmax) && (imp > ceps)){
      it=it+1  
      kk=0
      H = matrix(0,ndim,ndim)
      T0 = matrix(0,n,ndim)
      
      for(j in 1:q){
        k = kk+1
        kk = kk + zncati[j]
        Zk = z[,k:kk]
        Dz = diag(rr[k:kk])
        ZkFv = t(Zk) %*% Fv
        W = solve(Dz) %*% ((ZkFv))
  #      op = ordinalY(rr[k:kk],W,1,itermax=100,eps=1e-6,verbose=0)
  #      W = op$yhat
        H = H + t(W) %*% Dz %*% W
        T0 = T0 + Zk %*% W
        A[k:kk,]=W
      }
      
      H=H/q 
      svdRes=svd(H)
      d=diag(svdRes$d)
      v=svdRes$v
      inv_d_sq=diag(1/sqrt(svdRes$d))
      inv_d_sq[which(inv_d_sq==Inf)] = 0
      Tv=(1/q) * T0 %*% v %*% inv_d_sq
      #A = A %*% v %*% inv_d_sq
      center = pseudoinverse(t(U) %*% U) %*% (t(U) %*% Tv)
      
      outK = kmeans(Tv,centers=center)
      center=outK$centers
      index = outK$cluster
      
      U=matrix(0,n,nclus)
      for(i in 1:n){
        U[i,index[i]]=1
      }
      
      Fv = U %*% center 
      svdRRes=svd (t(Fv) %*% Fv)
      dd=diag(svdRRes$d)
      inv_dd_sq = diag(1/sqrt(svdRRes$d))
      vv = svdRRes$v
      Fv = Fv %*% vv %*% inv_dd_sq
      #A = A%*% vv %*% diag(svdRRes$d)
      f=0
      kk=0
      for(j in 1:q){
        k=kk+1
        kk = kk + zncati[j]
        dif = Fv-z[,k:kk] %*% A[k:kk,]
        f=f + sum(diag(t(dif)%*% dif))
      }
      imp= f0 - f
      f0=f
    }# endwhile
    if (f <= oldf){
      oldf = f
      Uold = U
      Aold = A
      Fvold=Fv
      indexold = index
      centeroid = center
      
    }
  }
  f = oldf
  U = Uold
  A = Aold
  Fv=Fvold
  index = indexold
  center = centeroid
  group_membership = index
  nclass = apply(U,2,sum)
  out=list()
  out$obscoord=Tv
  out$attcoord=A
  out$centroid=center
  out$U = U
  out$cluID=index
  out$criterion=f  
  out
}