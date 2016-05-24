list.historylabels = function(t, t.max=15){

  if(t>t.max) stop(paste("argument t cannot be greater than t.max=",t.max,sep="")) 

  ph.list=vector(mode="list",length=0)
  
  for(j in 1:(t-1)){
    fac=vector(mode="list",length=j)
    lapply(fac[1:j], function(x) factor(c(0,1))) -> fac[1:j]
    matindex=expand.grid(fac)
    matindex=as.data.frame(t(matindex))
    lapply(matindex[1:ncol(matindex)],function(x) factor(x,levels=c("0","1"))) -> matindex[1:ncol(matindex)]
    ph.list=c(ph.list,as.list(matindex))
  }
  
  ph.list=c("",ph.list)
  
  out=lapply(ph.list,paste,collapse="")
  
  return(out)
}
