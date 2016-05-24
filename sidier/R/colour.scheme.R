colour.scheme<-function(def=NA, N=NA, colors=c("green2","red","yellow","blue","DarkOrchid1","gray51","chocolate","cyan4","saddle brown","aquamarine","chartreuse","chocolate1","DarkOrchid3","gray18","gold","DarkOrchid4","green4","gray29", "sienna3","tan1","blue4","limegreen","gray73","bisque3","deeppink","red4","OliveDrab4","gray95", "salmon","DeepPink4","green yellow","gray4","hot pink","pink2","dark orange","gold3"))
{
	if(is.na(def[1])|length(def)!=N)
		if(N<length(colors))
		{
		coloresAux<-colors
		colo<-coloresAux[1:N]
		}
		else
		{
		coloresAux<-colors()[sample(c(1,23,25:152,203:259,361:657),N)]
		colo<-coloresAux[1:N]
		}

	if(is.na(def[1])==FALSE & length(def)==N)
	colo<-def
colo
}

## HACERLE TB UNA SALIDA PARA COLORES DE FUENTE Y QUE SE PONGA EN BLANCO LOS QUE NO SE VEAN BIEN EN NEGRO!
