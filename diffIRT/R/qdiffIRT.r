qdiffIRT=function(p,ai,vi,ter,sd2A,sd2V,maxRT,A,W,model,eps){
  ou=c()
  for(ii in 1:length(p)){
    tmp=function(x) pdiffIRT(x,ai,vi,ter,sd2A,sd2V,A,W,model,eps) -p[ii]
    if(ii==1) ou[ii] = uniroot(tmp,c(exp(ter)+0.001,maxRT))$root
      else ou[ii] = uniroot(tmp,c(ou[ii-1],maxRT))$root
    }
  return(ou)
}

