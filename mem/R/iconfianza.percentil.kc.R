iconfianza.percentil.kc <-
function(datos,q=0.50,nivel=0.95,ic=T,colas=2){
	x<-datos[!is.na(datos)]
	n<-length(x)
	if (n!=0){
  	pnor<-qnorm((1-nivel)/colas,lower.tail=FALSE)
  	j<-max(1,floor(n*q-pnor*sqrt(n*q*(1-q))))
  	k<-min(n,ceiling(n*q+pnor*sqrt(n*q*(1-q))))
  	if (floor(n*q)==n*q){
  		p<-(sort(x)[n*q] +sort(x)[n*q+1] )/2
  	}else{
  		p<-sort(x)[1+floor(n*q)] 
  	}
  	med<-p
  	l.i<-sort(x)[j]
  	l.s<-sort(x)[k]
	}else{
    med<-NA
    l.i<-NA
    l.s<-NA
  }
	if (ic) return(c(l.i,med,l.s)) else return(rep(med,3))
}
