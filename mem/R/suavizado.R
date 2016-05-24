suavizado <-
function(i.datos,hsuav=-1){
	m<-length(i.datos)
	datosx<-(1:m)[!is.na(i.datos)]
	datosy<-i.datos[!is.na(i.datos)]
	if (hsuav==-1) hsuav<-h.select(datosx,datosy)
	smr<-sm.regression(datosx,datosy,h=hsuav,display="none",eval.points=1:m)
	salida<-smr$estimate
	salida[salida<0]<-0
	salida[is.nan(salida)]<-NA
	return(salida)
}
