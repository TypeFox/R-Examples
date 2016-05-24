influenza <-
function(i.datos,i.tipo=2,i.nivel=0.95,i.tipo.curva=2,i.nivel.curva=0.95,i.tipo.umbral=5,i.nivel.umbral=0.95,i.n.max=-1,i.colas=1,i.tipo.boot="normal",i.iteraciones.boot=10000,i.metodo=2,i.parametro=2.8,i.niveles=c(0.40,0.90,0.975),i.n.max.temp=10){

  #i.datos<-cyl
  #i.tipo=2
  #i.nivel=0.95
  #i.tipo.curva=2
  #i.nivel.curva=0.95
  #i.tipo.umbral=2
  #i.nivel.umbral=0.95
  #i.n.max=0
  #i.colas=1
  #i.tipo.boot=NA
  #i.iteraciones.boot=NA
  #i.metodo=1
  #i.niveles=c(0.40,c.90,0.975)
  #i.n.max.temp=-1
  
  if (is.matrix(i.datos)) i.datos<-as.data.frame(i.datos)
  
  if (i.n.max.temp>0) i.datos<-i.datos[(max((dim(i.datos)[2])-i.n.max.temp+1,1)):(dim(i.datos)[2])]
  
	datos<-apply(i.datos,2,fill.missing)
	
	semanas<-dim(datos)[1]
	anios<-dim(datos)[2]
	
  # Calcular el optimo

  if (i.n.max==-1){
    n.max<-max(1,round(30/anios,0))
    #if (anios>=10) n.max=3 else n.max=5
  }else if (i.n.max==0){
    n.max=semanas
  }else{
    n.max=i.n.max
  }
  
  optimo<-apply(datos,2,localizar.epidemia,i.n.values=n.max,i.metodo=i.metodo,i.parametro=i.parametro)

  datos.duracion.real<-extraer.datos.optimo.map(optimo)

	# por defecto la aritmetica en ambos
	ic.duracion<-iconfianza(as.numeric(datos.duracion.real[1,]),i.nivel,i.tipo,T,i.tipo.boot,i.iteraciones.boot,2)
  ic.duracion<-rbind(ic.duracion,c(floor(ic.duracion[1]),round(ic.duracion[2]),ceiling(ic.duracion[3])))	
	ic.porcentaje<-iconfianza(as.numeric(datos.duracion.real[2,]),i.nivel,i.tipo,T,i.tipo.boot,i.iteraciones.boot,2)

	duracion.media<-ic.duracion[2,2]
	
	########################################################################################
	# Calcula todos los parametros asociados a la temporada gripal en base a la estimacion #
	# del numero de semanas de duracion total.                                             #
	########################################################################################

	# Datos de la gripe en la temporada REAL y en la temporada MEDIA, esta ultima sirve para unir las temporadas.
	
	semana.inicio<-as.integer(rownames(i.datos)[1])
	gripe<-array(dim=c(7,anios,2))
	
	datos.duracion.media<-extraer.datos.curva.map(optimo,duracion.media)
	
	for (j in 1:anios){
		## Semana en la q comienza la temporada.real
		gripe[1,j,1]<-datos.duracion.real[4,j]
		## Semana en la q acaba
		gripe[2,j,1]<-datos.duracion.real[5,j]
		## numero de casos en la temporada.real de gripe
		gripe[3,j,1]<-datos.duracion.real[3,j]
		## Numero de casos totales de la temporada global
		gripe[4,j,1]<-sum(datos[,j],na.rm=TRUE)
		## Porcentaje cubierto de casos de la temporada gripal
		gripe[5,j,1]<-datos.duracion.real[2,j]
		## Casos Acumulados justo antes de comenzar la temporada gripal
		if (gripe[1,j,1]>1) gripe[6,j,1]<-sum(datos[1:(gripe[1,j,1]-1),j],na.rm=TRUE) else gripe[6,j,1]<-0
		## Casos Despues de la temporada gripal
		if (gripe[2,j,1]<semanas) gripe[7,j,1]<-sum(datos[(gripe[2,j,1]+1):semanas,j],na.rm=TRUE) else gripe[7,j,1]<-0
		
		## Semana en la q comienza la temporada.media
		gripe[1,j,2]<-datos.duracion.media[4,j]
		## Semana en la q acaba
		gripe[2,j,2]<-datos.duracion.media[5,j]
		## numero de casos en la temporada.media de gripe
		gripe[3,j,2]<-datos.duracion.media[3,j]
		## Numero de casos totales de la temporada global
		gripe[4,j,2]<-sum(datos[,j],na.rm=TRUE)
		## Porcentaje cubierto de casos de la temporada gripal
		gripe[5,j,2]<-datos.duracion.media[2,j]
		## Casos Acumulados justo antes de comenzar la temporada gripal
		if (gripe[1,j,2]>1) gripe[6,j,2]<-sum(datos[1:(gripe[1,j,2]-1),j],na.rm=TRUE) else gripe[6,j,2]<-0
		## Casos Despues de la temporada gripal
		if (gripe[2,j,2]<semanas) gripe[7,j,2]<-sum(datos[(gripe[2,j,2]+1):semanas,j],na.rm=TRUE) else gripe[7,j,2]<-0
	}
	
	# Indice de estado temporal. Si es un 1 me dice que estamos en epoca pretemporada,
	# si pone un 2 indica que estamos en temporada gripal. Si pone un 3 indica que
	# estamos en posttemporada.

	indices.temporada<-array(dim=c(semanas,anios,2))
	
	for (j in 1:anios){
    indices.temporada[-((datos.duracion.real[4,j]):semanas),j,1]<-1
    indices.temporada[(datos.duracion.real[4,j]):(datos.duracion.real[5,j]),j,1]<-2
    indices.temporada[-(1:(datos.duracion.real[5,j])),j,1]<-3
    indices.temporada[-((datos.duracion.media[4,j]):semanas),j,2]<-1
    indices.temporada[(datos.duracion.media[4,j]):(datos.duracion.media[5,j]),j,2]<-2
    indices.temporada[-(1:(datos.duracion.media[5,j])),j,2]<-3
  }
    	
	## Estimacion del numero de semana en el  comienza la gripe
		
	ic.inicio<-iconfianza(as.numeric(datos.duracion.real[4,]),i.nivel,i.tipo,T,i.tipo.boot,i.iteraciones.boot,2)
  ic.inicio<-rbind(ic.inicio,c(floor(ic.inicio[1]),round(ic.inicio[2]),ceiling(ic.inicio[3])))
  ic.inicio<-rbind(ic.inicio,c(semana.absoluta(ic.inicio[2,1],semana.inicio),semana.absoluta(ic.inicio[2,2],semana.inicio),semana.absoluta(ic.inicio[2,3],semana.inicio)))
  ic.inicio<-ic.inicio[c(2,3,1),]
	inicio.medio<-ic.inicio[1,2]
	
	# Grafico del esquema de las temporadas gripales para su union
	
	longitud.esquema<-semanas+max.fix.na(gripe[1,,2])-min.fix.na(gripe[1,,2])
  inicio.epidemia.esquema<-max.fix.na(gripe[1,,2])
  fin.epidemia.esquema<-max.fix.na(gripe[2,,2])
  
	esquema.temporadas<-array(dim=c(longitud.esquema,anios+2,3))
	
	for (j in 1:anios){
		for (i in 1:semanas){
			diferencia<-inicio.epidemia.esquema-gripe[1,j,2]
			esquema.temporadas[i+diferencia,j,1]<-i
			esquema.temporadas[i+diferencia,j,2]<-semana.absoluta(i,semana.inicio)
			esquema.temporadas[i+diferencia,j,3]<-i.datos[i,j]
		}
	}
	
	diferencia<-inicio.epidemia.esquema-inicio.medio
	for (i in 1:semanas){
		if ((i+diferencia)>=1 & (i+diferencia)<=longitud.esquema){
			esquema.temporadas[i+diferencia,anios+1,1]<-i
			esquema.temporadas[i+diferencia,anios+1,2]<-semana.absoluta(i,semana.inicio)
		}
	}

  esquema.temporadas[,anios+2,c(1,2)]<-3
  esquema.temporadas[1:(inicio.epidemia.esquema-1),anios+2,1]<-1
 	esquema.temporadas[(inicio.epidemia.esquema):(inicio.epidemia.esquema+duracion.media-1),anios+2,1]<-2
 	esquema.temporadas[1:(inicio.epidemia.esquema-1),anios+2,2]<-1
 	esquema.temporadas[(inicio.epidemia.esquema):(inicio.epidemia.esquema+duracion.media-1),anios+2,2]<-2
	
	## Temporadas moviles
	
	temporadas.moviles<-esquema.temporadas[,c(1:anios),3]
			
	## Limites de la temporada
	
	limites.temporada<-c(inicio.medio,inicio.medio+duracion.media-1)
	limites.temporada<-rbind(limites.temporada,c(semana.absoluta(inicio.medio,semana.inicio),semana.absoluta(inicio.medio+duracion.media-1,semana.inicio)))
	
	# Limites relativos

	limites.esquema<-c(inicio.epidemia.esquema,inicio.epidemia.esquema+duracion.media-1)
	
	## En ttotal reordenamos la gripe para q todas las temporadas comiencen la semana "estinicio" (es decir, quitamos informacion
	## por delante y por detras de esquema.temporada[,,3] para que queden 35 semanas (las que hay).
	
	temporadas.moviles.recortada<-array(dim=c(semanas,anios))
	
	## Temporada gripal comienza "estinicio" y termina en "estinicio+temporada-1"
	
	diferencia<-inicio.medio-limites.esquema[1]
	
	for (j in 1:anios){
		for (i in 1:semanas){
			if ((i-diferencia)>=1 & (i-diferencia)<=longitud.esquema){
				temporadas.moviles.recortada[i,j]<-temporadas.moviles[i-diferencia,j]
			}else{
				temporadas.moviles.recortada[i,j]<-NA
			}
		}
	}
	
	## Dibujamos la curva epidemica tipo

	# Importante: ¿Quitamos los ceros?
	iy<-(is.na(temporadas.moviles.recortada) | temporadas.moviles.recortada==0)
	temporadas.moviles.recortada.no.ceros<-temporadas.moviles.recortada
	temporadas.moviles.recortada.no.ceros[iy]<-NA
 	curva.tipo<-t(apply(temporadas.moviles.recortada.no.ceros,1,iconfianza,nivel=i.nivel.curva,tipo=i.tipo.curva,ic=T,tipo.boot=i.tipo.boot,iteraciones.boot=i.iteraciones.boot,colas=2))	

	## Seleccionamos los periodos pre, epidemia, y post
	
	## PRE y POST-TEMPORADA GRIPAL
	
	## Como no se registra mas q de la semana 40 a la 20, tenemos q muchas de las semanas tienen tasa 0, eliminamos esas semanas,
	## ya q podrian llevar a una infraestimacion de la tasa base fuera de temporada
	
	pre.post.datos<-rbind(as.vector(as.matrix(extraer.datos.pre.epi(optimo))),as.vector(as.matrix(extraer.datos.post.epi(optimo))))
  epi.datos<-as.vector(as.matrix(apply(datos,2,max.n.valores,n.max=n.max)))	 	
	# IC de la linea basica de pre y post temporada

	# por defecto estaba la geometrica
	pre.d<-pre.post.datos[1,!(is.na(pre.post.datos[1,]) | pre.post.datos[1,]==0)]
	pre.i<-iconfianza(pre.d,nivel=i.nivel.umbral,tipo=i.tipo.umbral,ic=T,tipo.boot=i.tipo.boot,iteraciones.boot=i.iteraciones.boot,colas=i.colas)
	post.d<-pre.post.datos[2,!(is.na(pre.post.datos[2,]) | pre.post.datos[2,]==0)]
	post.i<-iconfianza(post.d,nivel=i.nivel.umbral,tipo=i.tipo.umbral,ic=T,tipo.boot=i.tipo.boot,iteraciones.boot=i.iteraciones.boot,colas=i.colas)
	pre.post.intervalos<-rbind(pre.i,post.i)
	epi.intervalos<-numeric()
	#niveles<-c(0.50,0.90,0.95)
	for (niv in i.niveles){
 	  epi.intervalos<-rbind(epi.intervalos,c(niv,iconfianza(epi.datos,nivel=niv,tipo=6,ic=T,colas=i.colas)))
 	}
   	
	## Ademas, añadimos las estimaciones de las lineas basicas antes y despues
	
	#lineas.basicas<-array(dim=c(semanas,3))
	#
	#for (i in 1:(inicio.medio-1)){
	#	lineas.basicas[i,]<-pre.post.intervalos[1,1:3]
	#}
	#
	#for (i in ((inicio.medio+duracion.media):semanas)){
	#	lineas.basicas[i,]<-pre.post.intervalos[2,1:3]
	#}

	resultados<-list(
    i.datos=i.datos,
    datos=datos,
		parametro.tipo=i.tipo,
		parametro.nivel=i.nivel,
		parametro.tipo.curva=i.tipo.curva,
		parametro.nivel.curva=i.nivel.curva,
		parametro.tipo.umbral=i.tipo.umbral,
		parametro.nivel.umbral=i.nivel.umbral,
    parametro.n.max=i.n.max,
    parametro.colas=i.colas,
    parametro.tipo.boot=i.tipo.boot,
    parametro.iteraciones.boot=i.iteraciones.boot,
    parametro.metodo=i.metodo,
    parametro.parametro=i.parametro,
    parametro.niveles=i.niveles,
    parametro.n.max.temp=i.n.max.temp,
    semana.inicio=semana.inicio,
    n.temporadas=anios,
    n.semanas=semanas,
    n.max=n.max,
    optimo=optimo,
    datos.duracion.real=datos.duracion.real,
    ic.duracion=ic.duracion,
		duracion.media=duracion.media,
		ic.porcentaje=ic.porcentaje,
		ic.inicio=ic.inicio,
    inicio.medio=inicio.medio,    		
		gripe=gripe,
		indices.temporada=indices.temporada,
		esquema.temporadas=esquema.temporadas,
		temporadas.moviles=temporadas.moviles,
		temporadas.moviles.recortada=temporadas.moviles.recortada,
		curva.tipo=curva.tipo,
		pre.post.datos=pre.post.datos,
		epi.datos=epi.datos,
		pre.post.intervalos=pre.post.intervalos,
    epi.intervalos=epi.intervalos)
	return(resultados)
}
