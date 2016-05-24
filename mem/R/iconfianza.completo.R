iconfianza.completo <-
function(i.datos,tipo=1,i.nivel=0.95){
  datos<-as.numeric(i.datos)
  datos[datos==-Inf]<-NA
  datos[datos==Inf]<-NA
  n<-sum(!is.na(datos))
  resultados<-list(media=NA)
  if (n!=0){
    med<-mean(datos,na.rm=T)
    if (tipo==1) std<-sqrt(var(datos,na.rm=T)/n) else std<-sqrt(var(datos,na.rm=T))
    resultados$media<-med
    pnor.2<-qnorm((1-i.nivel)/2,lower.tail=F)
    pnor.1<-qnorm((1-i.nivel),lower.tail=F)
    resultados$dos.colas<-c(med-pnor.2*std,med+pnor.2*std)
    resultados$mayor<-c(med-pnor.1*std,Inf)
    resultados$menor<-c(-Inf,med+pnor.1*std)
    
  }
  return(resultados)
}
