isoph.td=function(START, STOP, STATUS, Z, X, shape, K, maxdec, maxiter, eps){

  #data set with observed failure time
  oft.index=(which(START==0)-1)      #observed failure time (oft) for each subj
  oft.index=c(oft.index[-1], length(START)) #-1 for the first subj, with combining oft of the last subj
  
  TIME.oft=STOP[oft.index]
  STATUS.oft=STATUS[oft.index]
  Z.oft=Z[oft.index]
  n=length(STATUS.oft)
  if( sum(STATUS.oft)!=sum(STATUS) )
    stop("It is not allowed to have an event (status=1) other than the last interval of each subject.")

  #sorted by z
  order.z=order(Z.oft)
  t=TIME.oft[order.z]
  
    #create ZT
    t.obs=sort(unique(t))
    nt=length(t.obs)    

    ZT=matrix(0,n,nt) #ZT=ZT=matrix(NA,n,nt) does not work for the c code

    for(i in 1:n)
      for(j in oft.index[i]:(max(0,oft.index[i-1])+1))
        ZT[i,which(t.obs<=STOP[j])]=Z[j]
  
  zt=ZT[order.z,]
  status=STATUS.oft[order.z]
  z=sort(Z.oft)
  z.obs=unique(z[status==1])
  m=length(z.obs)
  
  #anchor
  k=sum(z.obs<K)
  if(k==0) k=1 #for the case when min(z.obs) < K  
  zk=z.obs[k]
  
  #counting process
  Y=dN=matrix(0,n,nt)       #row is the subj, col is the time corresponding z_(i);  
  
  for(i in 1:n){
    rank.t=which(t[i]==t.obs)
    Y[i,][1:rank.t]=1
    if(status[i]==1) dN[i,][rank.t]=1
  }  
  
  #RPA (with interval), Y, dN
  Y2=matrix(0,m,nt)
  dN2=matrix(0,m,nt)
  
  rpa.Y=.C('RPA_td', as.integer(n), as.integer(nt), as.integer(m), as.double(zt), as.double(z.obs) , as.double(Y), as.integer(dN), Y2=as.double(Y2), dN2=as.integer(dN2) )
  dN2=matrix(rpa.Y$dN2,m,nt)
  Y2 =matrix(rpa.Y$Y2, m,nt)
    
  #initial value
  try(beta.hat<-coxph(Surv(START,STOP,STATUS)~Z)$coefficient, silent=T)
  if(!is.numeric(beta.hat)) beta.hat=0.01

  if(shape=='increasing'){;       psi= abs(beta.hat)*(z.obs-zk)
  }else if(shape=='decreasing'){; psi=-abs(beta.hat)*(z.obs-zk);  }  
  
  #picm
  dNsum=colSums(dN2)
  Delta=rowSums(dN2)     

  #interval (RPA)
  int=list()
  int[[1]]=c(-Inf,z.obs[2])
  int[[m]]=c(z.obs[m],Inf)
  for(i in 2:(m-1))
    int[[i]]=c(z.obs[i],z.obs[i+1])
  
  #picm & beta for newton raphson algo
  iter=0;  dist=1;  beta=0

  if(is.null(X)){ #no trt group  
    dist=0; exp.beta=NA
    picm=picm.ft(psi,m,z.obs,zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
    if(picm$conv==0) stop
    psi.new=picm$psi.new
  }else{
    iter=0;  dist=1;  beta=0  
    
    X.oft=X[oft.index]
    x=X.oft[order.z]
    while(dist>=eps){  
      iter=iter+1
      if(iter>maxiter) break   
    
      #estimate psi
      picm=picm.ft(psi,m,z.obs,zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
      if(picm$conv==0) stop
      psi.new=picm$psi.new
      psi.full=BTFft(m, n, int, z, psi.new)    
      
      #estimate beta (Y1&x1 or Y2&x2 should be the same);
      beta.new=NR.ft(x,beta,psi.full,n,nt,Y,dN)
      
      #update;
      dist=sqrt(sum(psi.new-psi)^2)+(beta.new-beta)^2
      
      psi=psi.new
      beta=beta.new
    }
    exp.beta=round(exp(beta.new), maxdec)
    exp.beta=formatC( exp.beta, format='f', digits=maxdec)    
  }
  
  #picm result
  conv="converged"
  if(dist>=eps) conv="not converged"  

  
  #back to full rank (later)
  psi.obs=round(psi.new, maxdec)
  
  #level sets
  psi.uniq=unique(psi.obs)
  n.lv=length(psi.uniq)
  
  lv.sets=c()
  zmin=formatC( round(min(z),maxdec), format='f', digits=maxdec)
  zmax=formatC( round(max(z),maxdec), format='f', digits=maxdec)
  
  if(n.lv==1){ #only one level sets
    lv.sets[1]=paste('[',zmin,',',zmax,']', sep='')
  }else{
    lv=c()
    for(i in 1:n.lv)
      lv[[i]]=formatC( round(z.obs[which(psi.obs==psi.uniq[i])],maxdec)[1], format='f', digits=maxdec)
    for(i in 1:n.lv){
      if(i==1){
        lv.sets[1]=paste('[',zmin,', ',lv[2],')', sep='')
      }else if(i==n.lv){
        lv.sets[i]=paste('[',lv[i][1],', ',zmax,']', sep='')
      }else{
        lv.sets[i]=paste('[',lv[i],', ',lv[i+1],')', sep='')
      }
    }
  }
  psi.hat=formatC( unique(psi.obs), format='f', digits=maxdec)
  HR.hat=formatC( unique(exp(psi.obs)), format='f', digits=maxdec)
  
  est=data.frame(psi.hat=psi.hat, HR.hat=HR.hat, lv.set=lv.sets)
  names(est)=c("psi.hat","exp(psi.hat)","level set of psi.hat")
  
  z.range=range(Z)
  
  return(list(est=est, exp.beta=exp.beta, conv=conv, psi=psi.obs, z=z.obs, z.range=z.range, K=K, shape=shape, n=n, nevent=sum(STATUS), njump=m))
}