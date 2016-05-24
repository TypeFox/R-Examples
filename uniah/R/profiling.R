profiling.ft=function(m,H,q,psi.list,conv.list){
  
  #profiling algo excluding 1% 
  lf=c()
  s=round(m/2*0.01,0) #spiking pt +/ 1%
  
  for(k in 1:m){
    if(conv.list[k]==1){
      j=(k-s):(k+s)
      j=j[j>0]
      j=j[j<=m]
      
      psi=psi.list[[k]]
      psi=psi-psi[1] #compared psis at k=1, since spiking may occur at mode
      psi=psi[-j]
      
      lf[k]=t(psi)%*%H[-j,-j]%*%psi/2-t(psi)%*%q[-j]
    }
  }
  
  mode=which.min(lf) #na.rm is not necessary
  
  return(list(mode=mode,lf=lf))
}