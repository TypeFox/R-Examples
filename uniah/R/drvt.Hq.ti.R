drvt.Hq.ti=function(t.obs,dN,Y,n,m,nt, z,z.obs,x){
  #1. H & q
  dt=diff(c(0,t.obs))
  dNsum=colSums(dN)
  Ysum=colSums(Y)  

  H=matrix(0,n,n)
  q=matrix(0,n)
  
  diag(H)=colSums(t(Y)*dt) #or rowSums(t(t(Y)*dt))
  q=rowSums(dN)
  
  diag(H)=diag(H)-rowSums(t((t(Y)/Ysum)*dt))
  q=q-rowSums(t((t(Y)*dNsum)/Ysum))    

  #time consuming, but there is no way to handle it. C code is also slow.
  for(i in 1:(n-1))
    for(j in (i+1):n)
      H[i,j]=H[j,i]=-sum(Y[i,]*Y[j,]/Ysum*dt)
  
  #RPA
  H2=matrix(0,m,m)
  q2=matrix(0,m)
  
  int=list()
  int[[1]]=which(z<z.obs[2])
  int[[m]]=which(z.obs[m]<=z)
  for(i in 2:(m-1))
    int[[i]]=which(z.obs[i]<=z & z<z.obs[i+1])  
  
  for(i in 1:m){
    s=int[[i]]
    H2[i,i]=sum(H[s,s])
    q2[i]=sum(q[s])
  }
  
  for(i in 1:(m-1))
    for(j in (i+1):m){
      s=int[[i]]
      t=int[[j]]
      H2[i,j]=H2[j,i]=sum(H[s,t])
    }

  if(is.null(x)){
    return(list(H=H2, q=q2))
  }else{
    #2. H.circ, q.circ (p=1)
    x.bar=c()
    for(j in 1:nt)
      x.bar[j]=sum(Y[,j]*x[j])/Ysum[j]
    
    H.circ=rep(0,nt)
    q.circ=0
    for(i in 1:n){ #at each time point as a vector
      H.circ=H.circ+Y[i,]*(x[i]-x.bar)^2*dt
      q.circ=q.circ+(x[i]-x.bar)*dN[i,]
    }
    H.circ=sum(H.circ)
    q.circ=sum(q.circ)
    
    #3. H.diam (p=1)
    H.diam=rep(0,n)
    for(i in 1:n)
      H.diam[i]=sum(H.diam[i]+Y[i,]*(x[i]-x.bar)*dt)  
    
    #RPA for H.diam
    H.diam2=rep(0,m)  
    
    for(i in 1:m){
      s=int[[i]]
      H.diam2[i]=sum(H.diam[s])
    }  
    return(list(H=H2, q=q2, H.circ=H.circ, q.circ=q.circ, H.diam=H.diam2))
  }
}


     

      
