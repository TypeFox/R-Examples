fill.missing <-
function(i.datos){
	idatos.smo<-suavizado(i.datos,-1)
	donde.missing<-missings.inside(i.datos)
	idatos.fin<-i.datos
	idatos.fin[donde.missing]<-idatos.smo[donde.missing]
	return(idatos.fin)
}
