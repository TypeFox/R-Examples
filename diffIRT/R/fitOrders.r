fitOrders=function(nit,order){
  ou=matrix(0,fitCumChoose(nit,order),nit+1)
  r=1
    for(i in 1:order){
    	tmp=t(combn(nit,i))
    	nr=nrow(tmp)
      ou[r:(r+nr-1),1:(i+1)]=cbind(i,tmp)
      r=r+nr
    }
  return(ou)
}

