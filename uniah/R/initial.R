initial.unimodal.ft=function(gamma,z.obs,k,m){
  if(k==1){;        psi= abs(gamma)*(z.obs-z.obs[1]) 
  }else if(k==m){;  psi=-abs(gamma)*(z.obs-z.obs[m])
  }else{;           psi.lt=-abs(gamma)*z.obs[1:k];
        psi.rt= abs(gamma)*z.obs[(k+1):m];  psi.lt=psi.lt-psi.lt[k];
        psi.rt=psi.rt-psi.rt[1]
        psi=c(psi.lt,psi.rt)
  }
  
  return(psi)
}
  
initial.ushape.ft=function(gamma,z.obs,k,m)  {
  if(k==1){;        psi=-abs(gamma)*(z.obs-z.obs[1])
  }else if(k==m){;  psi= abs(gamma)*(z.obs-z.obs[m])
  }else{;           psi.lt= abs(gamma)*z.obs[1:k];
        psi.rt=-abs(gamma)*z.obs[(k+1):m];  psi.lt=psi.lt-psi.lt[k];
        psi.rt=psi.rt-psi.rt[1]
        psi=c(psi.lt,psi.rt)
  }
  return(psi)
}