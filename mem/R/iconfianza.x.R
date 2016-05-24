iconfianza.x <-
function(datos,nivel=0.95,ic=T,colas=2){
	datos[datos==-Inf]<-NA
	datos[datos==Inf]<-NA
	n<-sum(!is.na(datos))
	if (n!=0){
		pnor<-qnorm((1-nivel)/colas,lower.tail=FALSE)
		med<-mean(datos,na.rm=T)
		des<-sqrt(var(datos,na.rm=T))
		l.i<-med-pnor*des
		l.s<-med+pnor*des
	}else{
    med<-NA
    l.i<-NA
    l.s<-NA
  }
	if (ic) return(c(l.i,med,l.s)) else return(rep(med,3))
}
