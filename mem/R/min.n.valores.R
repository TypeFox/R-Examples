min.n.valores <-
function(mis.datos,n.min=1){
	mis.datos[mis.datos==Inf]<-NA
	mis.datos[mis.datos==-Inf]<-NA
	ordenado<-sort(mis.datos,decreasing=F)
	if (n.min>0) resultado<-ordenado[1:n.min] else resultado<-ordenado
	return(resultado)
}
