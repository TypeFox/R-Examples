#newton-raphson algo with p=1
NR.ft=function(x,beta,psi.full,n,nt,Y,dN, maxiter,eps){
  
  iter=0
  dist=1
  while(dist>=eps){   
    iter=iter+1
    if(iter>maxiter) break    
    
    xb=x*beta  
    Yest=matrix(,n,n)
    S1=rep(0,n)    
    for(j in 1:nt){
      Yest[,j]=Y[,j]*exp(psi.full+xb)
      S1[j]=sum(Yest[,j]*x)
    }
    S2=S1^2
    S0=colSums(Yest)
    
    #0/0=0
    idx=which(S0>0)
    S2=S2[idx]
    S1=S1[idx]
    S0=S0[idx]
    dN=dN[,idx]
    
    #update beta;
    U=0;    I=0
    for(i in 1:n){
      U=U+sum( ((x[i]-S1/S0)*dN[i,]) ) #sum over time
      I=I+sum( ((-S2/S0+1)*dN[i,]) )
    }
    
    beta.new = beta - U/I    
    
    #distance
    dist=(beta.new-beta)^2
    beta=beta.new    
  }
  
  return(beta.new)
}