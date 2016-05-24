iconfianza.logx <-
function(datos,nivel=0.95,ic=T,colas=2){
	if (any(datos>0,na.rm=T)){
		lx<-log(datos)
		med<-exp(iconfianza.x(lx,nivel,ic,colas))
		return(med)
	}else{
		return(rep(NA,3))
	}
}
