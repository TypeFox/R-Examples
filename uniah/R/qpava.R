#1 qpava for unimodal
qpava.unimodal.ft=function(psi,z.obs,H2,h,q2,m,k,eps,maxiter){
  iter=0
  dist=1
  while(dist>eps){
    iter=iter+1    
    if(iter>maxiter) break
    
    #qpava
    g=(q2-H2%*%psi)/h
    if(sum(is.na(g))+sum(is.infinite(g))>=1){
      psi.new=NA
      break
    }

    #unimodal
    if(k==1){
      psi.new=-pava(-g,h)    
    }else if(k==m){
      psi.new= pava( g,h)        
    }else{
      psi.lt= pava( g[1:k],h[1:k]);
      psi.rt=-pava(-g[(k+1):m],h[(k+1):m])
      psi.new=c(psi.lt,psi.rt)
    }
  dist=sum(abs(psi-psi.new))
  psi=psi.new
  }
  #impose the anchor
  psi.new=psi.new-psi.new[k]  
  
  conv=0  
  if(dist<eps) conv=1
  
  return(list(psi.new=psi.new, conv=conv))
}

#2 qpava for ushape
qpava.ushape.ft=function(psi,z.obs,H2,h,q2,m,k,eps,maxiter){
  iter=0
  dist=1
  while(dist>eps){
    iter=iter+1    
    if(iter>maxiter) break
    
    #qpava
    g=(q2-H2%*%psi)/h
    if(sum(is.na(g))+sum(is.infinite(g))>=1){
      psi.new=NA
      break
    }
    
    #ushape
    if(k==1){
      psi.new= pava( g,h)    
    }else if(k==m){
      psi.new=-pava(-g,h)        
    }else{
      psi.lt=-pava(-g[1:k],h[1:k])
      psi.rt= pava( g[(k+1):m],h[(k+1):m])
      psi.new=c(psi.lt,psi.rt)
    }
    psi.new=psi.new-psi.new[k]
    
    dist=sum(abs(psi-psi.new))
    psi=psi.new
  }
  #impose the anchor
  psi.new=psi.new-psi.new[k]  
  
  conv=0
  if(dist<eps) conv=1
  
  return(list(psi.new=psi.new, conv=conv))
}

