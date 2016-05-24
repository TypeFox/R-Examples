###############################################
# rgf.date: 
# 
# 
# 
###############################################



###############################################
# NAME: rgf.summarize
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################
rgf.summarize=function(x){
	o=list(denx='',deny='',box='',summ='',vNA='',size=length(x))
	deny=denx=box=summ=0
	cont=1
	pb=txtProgressBar(min=0,max=o$size,char='*',width=20,style=3)
	for (i in x) {	
		b=readGDAL(i,silent=T)$band1
		d=density(b,n=75,na.rm=T)
		denx=cbind(denx,d$x)
		deny=cbind(deny,d$y)
		box=cbind(box,boxplot(b,plot=F)$stats)
		aux=summary(b)
		if (length(aux)!=7) aux=c(aux,"NA's"=0)
		summ=(cbind(summ,aux))	
		setTxtProgressBar(pb, cont)
		cont=cont+1
	}
	o$denx=denx[,-1]
	o$deny=deny[,-1]
	o$box=box[,-1]
	o$summ=summ[,-1]	
	o$vNA=o$summ[7,]
	o	
}

###############################################
# NAME: rgf.plot
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################
rgf.plot=function(o,type='rain'){
	if (type=='rain') {
		barplot(o$summ[6,],col=4,names.arg=1:200)
		barplot(o$summ[4,],add=T,col=3)
		barplot((o$summ[6,]==0)*-10,add=T,col='red')
	}
	if (type=='box') {
		boxplot(o$box)
	}
	if (type=='density') {
		plot(o$deny~o$denx,type='l')
	}	
}

###############################################
# NAME: rgf.create
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################

rgf.create = function(prefix,suffix='',ini,fin=ini,monthini=1,output) {
	if (missing(ini)) stop('ini parameter is missing')
	if (missing(fin)) {rgf.create(prefix,suffix,ini,ini,monthini,output)} 
	else {
		ini=as.numeric(ini)
		fin=as.numeric(fin)
		monthini=as.numeric(monthini)
		if (monthini>1) fin=fin+1 
		rg=paste(prefix,rep(ini:fin,each=9),'0',1:9,suffix,sep='')
		rg=c(rg,paste(prefix,rep(ini:fin,each=3),10:12,suffix,sep=''))
		rg=sort(rg)
		if (monthini>1) rg=rg[monthini:(length(rg)-(12-monthini+1))] 		
		if (!missing(output)) write(c(length(rg),rg),file=output)
		rg
	}
}

###############################################
# NAME: 
# PURPOSE: read a rgf file and convert it to a character vector
# INPUTS: file name
# OUTPUTS: a modified character vector of text in inFl,
#          1ª modificacion: si inFl 
###############################################

rgf.read = function (inFl){
	if (!file.exists(inFl)) stop(sprintf('file not found %s',inFl))
	if (file.info(inFl)['isdir'] == TRUE) stop(sprintf('file not found %s',inFl))
	dir=dirname(inFl)
	if (dir=='.') dir=getwd()
	er=try((aux=as.vector(read.table(inFl,sep='&')[,1])),silent=TRUE)
	if (class(er)=='try-error') stop (sprintf('error reading %s',inFl))
	if (!is.na(as.integer(aux[1]))) aux=aux[-1] #comprobar si hay un numero en la primera linea y eliminarlo
	if (basename(aux[1]) == aux[1]) aux=paste(dir,'/',aux,sep='')
	aux
}

###############################################
# NAME: rgf.when
# PURPOSE:
#     Compara una imagen de referencia con una lista de imagenes.
#     Para cada pixel se localiza la primera ocurrencia del valor 
#     de referencia en la lista de imagenes.
#
# INPUTS:
#       inFl: File names list of grid
#       ref: reference image 
#       order: 'FIRST' default, registra la primera ocurrencia
#              'LAST' registra la ultima ocurrencia
#       silent: logical Flag; if TRUE, comments outputs are supressed
# OUTPUTS:
#       A SpatialGridDataFrame/SpatialPixelDataframe class. 
# CHANGES: 20/04/2010 - añadir comprobacion de parametros 
# CHANGES: 20/04/2010 - bug en comportamiento de silent
# CHANGES: 20/04/2010 - bug en comportamiento de order
###############################################

rgf.when = function (inFl,ref,order='FIRST',silent=FALSE) {
	#comprobacion de parametros
	if (!(order %in% c('FIRST','LAST'))) stop('order must be FIRST or LAST')
	
	if (!silent) print(paste('Searching for',order,'ocurrences'))	
	#leer imagen de referencia
	ref=readGDAL(ref,silent=TRUE)
	
	#iniciar variables 
	n=length(inFl)	
	when=0
	aux=0
	if (!silent) pb =txtProgressBar(min=1,max=n,char='*',width=20,style=3)
	if (order=='FIRST') {range=n:1} #comparar imagenes de atras hacia adelante, para quedarme con la mas reciente
	else {range=1:n} #comparar de delante hacia atras...
	
	#bucle que recorre inFl buscando concidencias con ref 
	#a cada paso asigna a los pixeles coincidentes su posicion en la lista 
	#realiza un overlay transparente entre when y msk
	#las ultimas coincidencias encontradas, machacan a las previas...
	for (i in range) {
		aux=readGDAL(inFl[i],silent=TRUE)
		msk=(aux$band1==ref$band1)
		when=(when*!msk) + (i*msk) #overlay transparente entre when y msk 
		if (!silent) { 
			if (order=='FIRST') {setTxtProgressBar(pb, n-i+1)}
			else {setTxtProgressBar(pb,i)}
		}
	}
	if (!silent) close(pb)
	aux$band1=when
	aux
}






###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
# CHANGES:	27/01/2010 	- Bug en el conteo del % realizado para operaciones no acumulativas
# CHANGES:	07/02/2013 	- Bug en el calculo de MAX y MIN
###############################################

rgf.summary = function(inFl,outFl,step=length(inFl),fun='SUM',silent=FALSE,...) {
	
	if (length(inFl)%%step!=0) stop(paste('elements in inFl must be a multiply of',step))
	if (length(outFl)!=floor(length(inFl)/step)) stop('elements in outFl is not correct')
	if (step==1) stop('step must be greater than 1')
	
	#Extraer el numero de elementos del array 
	n=length(inFl)  #num de ficheros de entrada
	nah=floor(n/step) #num de ficheros de salida
	
	print(paste('Procesing', n, 'files. Summary', fun, 'every',step))
	
	#branch para funciones acumulativas o no acumulativas (hay que leerlas por trozos)
	if (fun %in% c('SUM','MAX','MIN','MEAN')) {
		
		if (!silent) pb =txtProgressBar(min=0,max=n,char='*',width=20,style=3)
		
		aux=readGDAL(inFl[1],silent=TRUE)	
		
		for(i in 0:(nah-1)){
			#CHANGE dat=0 
			#inicializo dat con 0 para SUM o MEAN o con la primera imagen del grupo para MAX o MIN
			switch (fun,
					SUM = {dat = 0},
					MIN = {dat = readGDAL(inFl[i*step+j],silent=TRUE)$band1},
					MAX = {dat = readGDAL(inFl[i*step+j],silent=TRUE)$band1},
					MEAN= {dat = 0},
			)
			for(j in 1:step) {
				switch (fun,
						SUM = {dat = dat+readGDAL(inFl[i*step+j],silent=TRUE)$band1},
						MIN = {dat = pmin(dat,readGDAL(inFl[i*step+j],silent=TRUE)$band1)},
						MAX = {dat = pmax(dat,readGDAL(inFl[i*step+j],silent=TRUE)$band1)},
						MEAN= {dat = dat+readGDAL(inFl[i*step+j],silent=TRUE)$band1/step},
				)
				if (!silent) setTxtProgressBar(pb, i*step+j)
			}
			aux$band1=dat
			writeGDAL(aux,outFl[i+1],...)
		}
	} else { #funciones no acumulativas 
		
		#num de filas de la imagen
		rows=GDALinfo(inFl[1])[1]
		#num de columnas de la imagen
		cols=GDALinfo(inFl[1])[2]
		
		#calculo size del buffer de lectura para que lea bloques de 100MB
		linesToRead=ceiling(100000000/(cols*n*8))
		#num de bloques de size LinesToRead en la imagen
		nblocks=ceiling(rows/linesToRead)
		
		if (!silent) pb =txtProgressBar(min=0,max=nah*nblocks-1,char='*',width=20,style=3)
		
		for(i in 0:(nah-1)){
			dat=0
			for (k in 0:(nblocks-1)){
				outdf=0
				for(j in 1:step) {
					aux=GDAL.open(inFl[i*step+j],TRUE)
					#offsets
					offsetIni=k*linesToRead+1
					offsetFin=offsetIni+linesToRead-1
					#no rebasar fin de archivo
					if (offsetFin>rows) offsetFin=rows
					#lee bloque de la imagen SIN cargarla completamente en memoria 
					Strip=aux[offsetIni:offsetFin,]
					#cerrar conexion
					GDAL.close(aux)
					#anexar como columna al dataframe de salida    
					outdf=cbind(outdf,Strip$band1)
				} 
				#actualizo progressbar		
				if (!silent) setTxtProgressBar(pb, i*nblocks+k)
				
				outdf=outdf[,-1]
				switch (fun,
						MEDIAN= {dat=c(dat,rowMedians(outdf))},
						RANGE = {dat=c(dat,rowRanges(outdf))},
						SD = {dat=c(dat,rowSds(outdf))},
						VAR= {dat=c(dat,rowVars(outdf))},
						COUNT ={dat=c(dat,rowCounts(outdf))}
				)					
			}
			
			aux=readGDAL(inFl[1],silent=TRUE)		
			dat=dat[-1]
			aux$band1=dat
			writeGDAL(aux,outFl[i+1],...)
		}
	}
	if (!silent) close(pb)
}
