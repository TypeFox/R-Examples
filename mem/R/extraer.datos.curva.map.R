extraer.datos.curva.map <-
function(i.epi,i.val){
  resultado<-as.data.frame(i.epi[[1]]$curva.map[i.epi[[1]]$curva.map[,1]==i.val])
  m<-ncol(resultado)
  n<-length(names(i.epi))
  if (n>1) for (i in 2:n){
    resultado<-cbind(resultado,i.epi[[i]]$curva.map[i.epi[[i]]$curva.map[,1]==i.val])
  }
  names(resultado)<-rep(names(i.epi),each=m)
  return(resultado)
}
