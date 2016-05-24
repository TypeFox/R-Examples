calcular.map <-
function(i.datos){
  datos<-as.vector(as.matrix(i.datos))
  semanas<-length(datos)
  maxsumasemanas<-array(dim=c(semanas,5))
  for (s in 1:semanas){
    sumasemanas<-numeric()
    for (i in 1:(semanas+1-s)){
      sumasemanas<-c(sumasemanas,sum(datos[i:(i+s-1)],na.rm=TRUE))
    }
    maxsumasemanas[s,1]<-s
    maxsumasemanas[s,3]<-max.fix.na(sumasemanas)
    maxsumasemanas[s,4]<-min((1:(semanas+1-s))[max.fix.na(sumasemanas)==sumasemanas])
    maxsumasemanas[s,5]<-maxsumasemanas[s,4]+s-1
  }
	sumaanual<-sum(datos,na.rm=TRUE)
	maxsumasemanas[,2]<-100*maxsumasemanas[,3]/sumaanual
	return(maxsumasemanas)
}
