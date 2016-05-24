uniah.ti.known=function(TIME, STATUS, Z, X, shape, K, maxdec, maxiter, eps){  
  null.x=0
  if(is.null(X)) null.x=1
  
  #sorted by z
  n=length(STATUS)
  order.z=order(Z)
  t=TIME[order.z]

  status=STATUS[order.z]
  z=sort(Z)
  z.obs=unique(z[status==1])
  m=length(z.obs)
  k=sum(z.obs<K)
  if(k==0) k=1
  
  t.obs=sort(unique(t))
  nt=length(t.obs)  
  
  #counting & at risk processes
  Y=dN=matrix(0,n,nt)       #row is the subj, col is the time corresponding z_(i);
  for(i in 1:n){
    rank.t=which(t[i]==t.obs)
    Y[i,][1:rank.t]=1
    if(status[i]==1) dN[i,][rank.t]=1
  }

  #Hessian with RPA
  if(null.x==1){  
    drvt.Hq=drvt.Hq.ti(t.obs,dN,Y,n,m,nt,z,z.obs,x=X)
  }else{
    x=X[order.z]
    drvt.Hq=drvt.Hq.ti(t.obs,dN,Y,n,m,nt,z,z.obs,x=x)
    
    H.circ=drvt.Hq$H.circ
    q.circ=drvt.Hq$q.circ
    H.diam=drvt.Hq$H.diam
  }
  H=drvt.Hq$H
  q=drvt.Hq$q
  
  #initial 
  t2=t+runif(n,0.01) #ad-hoc to avoid ties
  try(gamma<-summary(ahaz(Surv(t2,status), z))$coef[,1], silent=T) #for inital value
  if(!is.numeric(gamma)) gamma=0.01
  
  if(shape=="unimodal"){; psi=initial.unimodal.ft(gamma,z.obs,k,m)
  }else{;                 psi=initial.ushape.ft(gamma,z.obs,k,m);  }  
  
  #quadratic pava  
  h=diag(H)
  H2=H; diag(H2)=0 
  
  if(shape=="unimodal"){;  qpava.ft=qpava.unimodal.ft
  }else{;                  qpava.ft=qpava.ushape.ft;  }

  h=diag(H)
  H2=H; diag(H2)=0  

  if(null.x==1){ #no treatment effect
    beta=NA
    qpava=qpava.ft(psi,z.obs,H2,h,q, m,k,eps,maxiter)
    psi.new=qpava$psi.new
    conv="converged"
    if(qpava$conv==0) conv="not converged"    
  }else{
    conv="not converged"
    dist=1;  beta=0;  iter=0;    
    q2=q-c(H.diam)*beta
    
    while(dist>eps){
      iter=iter+1
      if(iter>maxiter) break
      
      #(a) psi.new
      qpava=qpava.ft(psi,z.obs,H2,h,q2, m,k,eps,maxiter)
      psi.new=qpava$psi.new
      conv=qpava$conv
      
      if(conv==0) stop("Algorithms were not converged.")
      
      #profiling algorithm with excluding 1% 
      
      #(b) beta.new (p=1)
      beta.new=c( 1/H.circ * (q.circ - t(H.diam)%*%psi.new) )
      
      #(c) update psi and beta
      dist=sum(abs(psi.new-psi))+abs(beta.new-beta)
      
      psi=psi.new
      beta=beta.new
    }
    
    if(iter>maxiter){; stop("Algorithms were not converged.")
    }else{; conv="converged";}
  }

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
  #beta=formatC( unique(beta), format='f', digits=maxdec)  

  est=data.frame(psi.hat=psi.hat, lv.set=lv.sets)
  names(est)=c("psi.hat", "level set")
  
  #for plot
  psi.obs=c(psi.obs[1],psi.obs,psi.obs[m])
  z.obs=c(min(z),z.obs,max(z))  

  #return(list(est=est, conv=conv, psi=psi.obs, z=z.obs, M=K, shape=shape, iter=iter, dist=dist, n=n, nevent=sum(STATUS), njump=m))
  return(list(est=est, beta=beta, conv=conv, psi=psi.obs, z=z.obs, M=K, shape=shape, n=n, nevent=sum(STATUS), njump=m))
}
