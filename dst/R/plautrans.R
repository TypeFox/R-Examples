plautrans<-function(x) {
  # Plausibility transformation on the distribution of singletons
  # 
  # x ix a list obtained from the tabresul function.
  # We extract the singletons which appears in the input table.
  # We compute the plausibility transformation on the singletons.
  # 
  zx<-x$mbp
  nsing<- -4+ncol(zx)
  z1<-apply(zx[,1:nsing],1,sum)
  zs<-matrix(zx[z1==1,],ncol=ncol(zx))
  # calculate distribution
  trplau<-zs[,nsing+3]/sum(zs[,nsing+3])
  y<-cbind(zs[,1:nsing],trplau)
  return(y)
}
