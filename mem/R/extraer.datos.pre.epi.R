extraer.datos.pre.epi <-
function(i.epi){
  resultado<-as.data.frame(i.epi[[1]]$pre.epi)
  m<-ncol(resultado)
  n<-length(names(i.epi))
  if (n>1) for (i in 2:n){
    resultado<-cbind(resultado,i.epi[[i]]$pre.epi)
  }
  names(resultado)<-rep(names(i.epi),each=m)
  return(resultado)
}
