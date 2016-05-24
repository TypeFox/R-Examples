
editr2dRfile = function(conf='') {
###############################################
# NAME: r2dRueWiz
# PURPOSE:
#     Read an 2dRue configuration file and perform actions acordely 
#     Can calculate the monitoring, the assessment, both or none, acordely with the input.
# INPUTS:
#       conf: config file. If none is provided, it will be create through interactive questions to the user 
# OUTPUTS:
#       r2dRue wiz is a procedure that can produced all the output involved in r2dRue analisys.
#		It will be created a log file named as the config file but with its extension turned to .log	
	input=function(tag,default=''){
		if (default!='') tag=paste(tag,' [',default,']',sep='')
		repeat {
			aux=readline(paste(tag,' ?: ',sep=''))
			if (aux!='') {
				return(aux);break
			} 
			if ((aux=='') & (default!='')) {
				return(default);break
			}			
		}		
	}
	
	year=function(x) {as.integer(format(x,'%Y'))}
	month=function(x) {as.integer(format(x,'%m'))}
	
	canReadRgf = function(x) {		
		if (!is.na(file.info(x)$isdir) & (file.info(x)$isdir==FALSE)) {			
			aux=rgf.read(x)
			if (!is.na(file.info(aux[1])$isdir) & (file.info(aux[1])$isdir==FALSE)){
				img=readGDAL(aux[1],silent=TRUE)			
				if (class(img)=='SpatialGridDataFrame')  {return(TRUE)}				
			}
		}
		return(FALSE)
	}
	
	############################## main
	response=list(comment='',viRgf='',rainRgf='',petRgf='',mHidro='',acum='',pOut='',sYear='',sMonth='',yIni='',yEnd='',driver='',flag='')
		
	text="##########################################\n############# r2dRue Wizard ##############\n#                                        #\n\n"
	
	#options(show.error.messages=FALSE)
	try({
		conaux=as.data.frame(readConfigFile(conf),stringsAsFactors=FALSE)
		response[names(conaux)]=conaux
	})
	if (response$pOut=='') response$pOut=getwd()
	if (response$driver=='') response$driver='RST'
	if (response$flag=='') response$flag=-999
	if (conf=='') conf=paste('rue',format(Sys.time(), "%d%b"),'.conf',sep='')
	
	cat(text)
	try({
		response$comment=input('Description of this run',response$comment)		
		repeat {
			aux=input('Output directory',response$pOut)
			sub('[/\\]$','',aux)
			if (!is.na(file.info(aux)$isdir) & (file.info(aux)$isdir!=FALSE)) {response$pOut=aux; break}
			else print('the directory does not exist...create first')
		}	
		repeat {
			aux=input('Vegetation Index raster group',response$viRgf)
			if (canReadRgf(aux)) {response$viRgf=aux; break}
			else print('cant read this raster group...')
		}
		repeat {
			aux=input('Precipitation raster group',response$rainRgf)			
			if (canReadRgf(aux)) {response$rainRgf=aux; break}
			else print('cant read this raster group...')
		}
		repeat {
			aux=input('PET raster group',response$petRgf)
			if (canReadRgf(aux)) {response$petRgf=aux; break}
			else print('cant read this raster group...')
		}
		repeat {
			aux=input('Start moment (yyyy/mm) of these raster groups',paste(response$sYear,response$sMonth,sep='/'))
			er=try({aux=as.Date(paste(aux,'01',sep='/'))})		
			if (class(er)!='try-error') {sDate=aux; break}
			else print('not a date')
		}
		repeat {
			aux=input('Start month of hydrological year [1-12]',response$mHidro)
			er=try({iaux=as.integer(aux)})		
			if (iaux %in% 1:12) {response$mHidro=iaux; break}
			else print('its simple... [1-12]')
		}
	})
	
	#options(show.error.messages=TRUE)
	svi=rgf.read(response$viRgf)
	srain=rgf.read(response$rainRgf)
	spet=rgf.read(response$petRgf)	
	
	#calculate dates of the elements in the serie 
	sLength=length(svi)
	sIniDate=sDate
	sDates=seq(sIniDate,length.out=sLength,by='month')
	sEndDate=sDates[sLength]
	
	sIniYear=year(sIniDate)
	sEndYear=year(sEndDate)
	
	sHYears=sFailHYears=0
	for (i in sIniYear:sEndYear) {
		aux1=as.Date(paste(i,response$mHidro,1,sep='/'))
		aux2=seq(aux1,length.out=12,by='month')[12]
		if ((aux1 %in% sDates) & (aux2 %in% sDates)) sHYears=c(sHYears,i)  
		else sFailHYears=c(sFailHYears,i)
	}	
	sHYears=sHYears[-1]
	sFailHYears=sFailHYears[-1]
		
	#show info
	cat(sprintf('\nOriginal data: %d images, from %s to %s',sLength,format(sIniDate,'%b/%Y'),format(sEndDate,'%b/%Y')))
	cat(sprintf('\n             : %d Hydrological years, %s - %s starting in %s',length(sHYears),sHYears[1],sHYears[length(sHYears)],month.name[response$mHidro]))
	cat(sprintf('\n             : %d Incomplete Hydrological years: %s\n',length(sFailHYears),paste(sFailHYears,collapse=', ')))
		
	repeat {
		aux=input('Number of cumulative months for preceding rain',response$acum)
		er=try({iaux=as.integer(aux)})		
		if (iaux %in% 1:30) {response$acum=iaux; break}
		else print('its simple... a integer great than 0')
	}
	
	preDate=sort(seq(as.Date(paste(sHYears[1],response$mHidro,1,sep='/')),length.out=response$acum+1,by='-1 months')[-1])
		
	if (!(preDate[1] %in% sDates)) sHYears=sHYears[-1]
	
	repeat {
		aux=input(sprintf('Start year of this run [%d-%d]',sHYears[1],sHYears[length(sHYears)]),response$yIni)
		er=try({iaux=as.integer(aux)})		
		if (iaux %in% sHYears) {response$yIni=iaux; break}
		else sprintf('[%d-%d]',sHYears[1],sHYears[length(sHYears)])
	}
	repeat {
		aux=input(sprintf('End year of this run [%d-%d]',sHYears[1],sHYears[length(sHYears)]),response$yEnd)
		er=try({iaux=as.integer(aux)})		
		if (iaux %in% sHYears) {response$yEnd=iaux; break}
		else sprintf('[%d-%d]',sHYears[1],sHYears[length(sHYears)])
	}
	repeat {		
		aux=input('GIS format for output images',response$driver)				
		if (isSupportedGDALFormat(aux)) {response$driver=aux; break}
		else print('a valid write GDAL driver')
	}
	repeat {
		aux=input('Missing value for output images',response$flag)
		er=try({iaux=as.integer(aux)})		
		if (is.integer(iaux)) {response$flag=iaux; break}
		else print('integer')
	}
	
	response$sYear=year(sIniDate)
	response$sMonth=month(sIniDate)
	
	fileName=input('File name for this config file',conf)	
	write.table(t(as.data.frame(response)),fileName,sep='=',quote=FALSE,col.names=FALSE)	
}

r2dRplot=function(o,type='rain',scope='run',var='vi',pixel=1,col=c('blue','green4','salmon','gray')){
##############################################
# NAME: ruePlot
# PURPOSE:
#     Read an 2dRue configuration file and perform actions acordely 
#     Caln calculate the monitoring, the assessment, both or none, acordely with the input.
# INPUTS:
	if (o$summarize) {
	switch(type,		
		'vimax'= {		
			nf <- layout(matrix(c(1,1,2,3,3,4), 3, 2))			
			a=readGDAL(o$viMax,silent=TRUE)
			image(a)
			title(main='vi maximum')
			plot(density(a$band1,na.rm=TRUE),main='',xlab='vi')
			a=readGDAL(o$viMaxWhen)
			image(a)
			title(main='when vi maximum')
			hist(o$sDates[a$band1],xlab='months',breaks='months',format='%b%y',las=3,freq=T)			
		},
		'box'= {
			if (var=='vi') boxplot(o$RESvi$box,names=format(o$rDates,'%b%y'),las=3,main='spatially lumped vegetation index',ylab='vegetation index')
			if (var=='rain') boxplot(o$RESrain$box,names=format(o$rDates,'%b%y'),las=3,main='spatially lumped precipitation',ylab='rain (mm)')		
		},
		'density'= {		
			if (var=='vi') plot(o$RESvi$deny~o$RESvi$denx,type='l',col=c(1:12))
			if (var=='rain') plot(o$RESrain$deny~o$RESrain$denx,type='l',col=c(1:12))
		},
		'rain'= {
			barplot(o$RESrain$summ[6,],col=4,names=format(o$rDates,'%b%y'),las=3,main='Precipitation',ylab='rain (mm)')
			barplot(o$RESrain$summ[4,],add=T,col=3)
			barplot((o$RESrain$summ[6,]==0)*-10,add=T,col='red')		
		},
		'assessment1'= {
			a=readGDAL(o$rueEx,silent=TRUE)
			a$rueObsEx=a$band1
			a$rueObsMe=readGDAL(o$rueMe,silent=TRUE)$band1
			a$aiObsEx=readGDAL(o$aiEx,silent=TRUE)$band1
			a$aiObsMe=readGDAL(o$aiMe,silent=TRUE)$band1
			nf <- layout(matrix(c(1,2), 1, 2))
			layout.show(nf)						
			plot(a$rueObsEx~a$aiObsEx,main='rue vs ai - Extremo',xlab='aiObsEx',ylab='rueObsEx')
			plot(a$rueObsMe~a$aiObsMe,main='rue vs ai - Medio',xlab='aiObsMe',ylab='rueObsMe')
			
		},
		'assessment2'= {
			a=readGDAL(o$rueEx,silent=TRUE)
			a$rueObsEx=a$band1
			a$rueObsMe=readGDAL(o$rueMe,silent=TRUE)$band1
			a$aiObsEx=readGDAL(o$aiEx,silent=TRUE)$band1
			a$aiObsMe=readGDAL(o$aiMe,silent=TRUE)$band1
			nf <- layout(matrix(c(1,2), 1, 2))			
			plot(spplot(a,zcol=c('rueObsEx','rueObsMe','aiObsEx','aiObsMe')))
		},
		'monitoring'= {
			a=readGDAL(o$f1,silent=TRUE)
			a$rueObsEx=a$f1band1
			a$rueObsMe=readGDAL(o$f2,silent=TRUE)$band1
			a$aiObsEx=readGDAL(o$f3,silent=TRUE)$band1
			a$f4=readGDAL(o$f4,silent=TRUE)$band1
			nf <- layout(matrix(c(1,2), 1, 2))
			layout.show(nf)
			plot(spplot(a,zcol=c('rueObsEx','rueObsMe','aiObsEx','aiObsMe')))
		},
		'pixel'={				
			pos=pixel*o$sLength*4
			in1=file(o$STKsvi,'rb')
			in2=file(o$STKsrain,'rb')
			seek(in1,pos)
			seek(in2,pos)
			bvi  =readBin(in1,numeric(),o$sLength,4)
			brain=readBin(in2,numeric(),o$sLength,4)		
			maxbrain=max(brain)			
			posmaxvi=which(bvi==max(bvi))
			#preposmaxvi=o$sDate[which(bvi==max(bvi))-o$acum]
			plot(o$sDate,bvi,type='l',ylim=c(0,1),col=col[2],ylab='vegetation index & rain (scaled)',xaxt='n',main=paste('VI & Rain for pixel ',pixel))
			axis.Date(1, at=seq(o$sIniDate,o$sEndDate,by='year'),las=3,cex.axis=1)
			# TODO: arreglar cuando - acum es negativo...
			if (length(posmaxvi)!=0) rect(o$sDates[posmaxvi-o$acum],0,o$sDates[posmaxvi],1,density=10,col=col[4],border=col[3])			
			lines(o$sDate,brain/maxbrain,type='h',col=col[1])
			abline(v=o$sDates[posmaxvi],col=col[3])
			lines(o$sDate,bvi,type='l',col=col[2],lwd=2)
			close(in1)
			close(in2)
		})
	} else (stop('summarize not done... ejecute summarize() function first'))
}

showInfo=function (o) {
###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################	
	#print info
	aux='\n'
	aux=c(aux,'\n',sprintf('################### r2dRue RUN: %s',o$comment))
	aux=c(aux,'\n',sprintf('Original data: %d images, from %s to %s',o$sLength,format(o$sIniDate,'%b/%Y'),format(o$sEndDate,'%b/%Y')))
	aux=c(aux,'\n',sprintf('Analysis data: %d images, from %s to %s, %d Hydrological years starting in %s',o$rLength,format(o$rIniDate,'%b/%Y'),format(o$rEndDate,'%b/%Y'),o$rYears,month.name[o$mHidro]))
	aux=c(aux,'\n',sprintf('               %d cumulative precipitation months, %d preceding images from  %s to %s',o$acum, length(o$ppet),format(o$rPreDates[1],'%b/%Y'),format(o$rPreDates[length(o$ppet)],'%b/%Y') ))
	aux=c(aux,'\n\n',sprintf('---------- Raster Group Info'))
	aux=c(aux,'\n',sprintf('viRgf     : %s',attr(o$vi,'path')))
	aux=c(aux,'\n',sprintf('            %s ...',paste(head(attr(o$vi,'files'),10),collapse=' ')))
	aux=c(aux,'\n',sprintf('rainRgf   : %s',attr(o$rain,'path')))
	aux=c(aux,'\n',sprintf('            %s ...',paste(head(attr(o$rain,'files'),10),collapse=' ')))
	aux=c(aux,'\n',sprintf('petRgf    : %s',attr(o$pet,'path')))
	aux=c(aux,'\n',sprintf('            %s ...',paste(head(attr(o$pet,'files'),10),collapse=' ')))
	aux=c(aux,'\n',sprintf('preRainRgf: %s',attr(o$prain,'path')))
	aux=c(aux,'\n',sprintf('            %s ...',paste(head(attr(o$prain,'files'),10),collapse=' ')))
	aux=c(aux,'\n',sprintf('prePetRgf : %s',attr(o$ppet,'path')))
	aux=c(aux,'\n',sprintf('            %s ...',paste(head(attr(o$ppet,'files'),10),collapse=' ')))
	aux=c(aux,'\n\n',sprintf('---------- Spatial Info'))
	aux=c(aux,'\n',sprintf('cols: %d  rows: %d  res: %f  proj: %s',o$gdal[1],o$gdal[2],o$gdal[6],attr(o$gdal,'projection')))
	aux=c(aux,'\n')
	cat(aux)
	if (o$assessment){
		cat(c('\n',paste('---------- Assessment results at ',attr(o$assessment,'date')),'\n'))		
		print(attr(o$assessment,'summary'))
	}else{cat(c('\n',sprintf('---------- assessment results not updated')))}
	if (o$monitoring){
		cat(c('\n',paste('---------- Monitoring results at ',attr(o$monitoring,'date')),'\n'))
		print(attr(o$monitoring,'summary'))
	}else{cat(c('\n',sprintf('---------- Monitoring results not updated')))}
	cat('\n')
	
}

###############################################
# NAME: 
# PURPOSE: Lee un fo
# INPUTS:
# OUTPUTS:
readConfigFile=function (conf) {	
	co=list(
			pOut='',mHidro='',acum='',viRgf='',rainRgf='',sYear='',sMonth='',
			yIni='',yEnd='',driver='',flag='',petRgf='',
			comment='',acction=''
	)
	obligatorios=names(co)[1:12]	
	f=readIniFile(conf)
	co[names(f)]=f	
	if (any(co[obligatorios] == '')) stop(sprintf('Faltan parametros obligatorios en %s, check the conf file',conf))	
	co
}

###############################################
# NAME: createRunConfig
# PURPOSE:
#     crea una lista con la informacion necesaria para una ejecucion de 2dRue
# INPUTS:
#     conf: lista de configuracion basica (sin campos calculados) 
# OUTPUTS:
readr2dRfile=function (conf){
	#campos forzados a enteros
	enteros=c('mHidro','acum','sYear','sMonth','yIni','yEnd')
	#definicion inicial
	o=list(
		#from config file
		comment='',pOut='',mHidro='',acum='',viRgf='',rainRgf='',sYear='',sMonth='',
		yIni='',yEnd='',driver='',flag='',acction='',petRgf='',
		#calculated		
		vi='',rain='',pet='',ppet='',prain='',rLength='',rIniDate='',rEndDate='',rDates='',rPreDates='',
		svi='',srain='',spet='',sIniDate='',sEndDate='',sDates='',sLength='',		
		#calculated spaciales
		gdal='',
		#calculates extern
		assessment=FALSE,rueMed='',rueEx='',aiMed='',aiEx='',
		monitoring=FALSE,f1='',f2='',f3='',f4='',f5='',f6='',f7='',f8='',f9='',
		summarize=FALSE,viMax='',whenviMax=''
	)
	
	#read fields from file
	co=readConfigFile(conf)
	o[names(co)]=co #copy fields from config file	
	o[enteros]=as.integer(co[enteros]) #pasar a enteros	
		
	#read rgf files
	o$svi=rgf.read(o$viRgf)
	o$srain=rgf.read(o$rainRgf)
	o$spet=rgf.read(o$petRgf)
	stopifnot(length(o$spet)==length(o$svi))
	
	#calculate dates of the elements in the serie 
	o$sLength=length(o$svi)
	o$sIniDate=as.Date(paste(o$sYear,o$sMonth,1,sep='/'))
	o$sDates=seq(o$sIniDate,length.out=o$sLength,by='month')
	o$sEndDate=o$sDates[o$sLength]
		
	o$rYears=o$yEnd-o$yIni+1
	o$rLength=o$rYears*12	
	o$rIniDate=as.Date(paste(o$yIni,o$mHidro,1,sep='/'))
	o$rDates=seq(o$rIniDate,length.out=o$rLength,by='month')
	o$rEndDate=o$rDates[o$rLength]
		
	o$rPreDates=sort(seq(o$rIniDate,length.out=o$acum+1,by='-1 month')[-1])

	if (!all(o$rDates %in% o$sDates)) {
		print(o$rDates)
		print(o$sDates)
		stop('check start and initial dates')
	}
	if (!all(o$rPreDates %in% o$sDates)) stop('check start and initial dates')
	
	#spatial atributes
	o$gdal=GDALinfo(o$svi[1],silent=TRUE)	
	if (!isSupportedGDALFormat(o$driver)) stop('not supported GDAL driver')
	if (!is.finite(as.numeric(o$flag))) stop('not a valid missing value flag ')
	
	#calculate vi,rain,pet,ppet and prain series
	o$vi=o$svi[o$sDates %in% o$rDates]
	o$rain=o$srain[o$sDates %in% o$rDates]
	o$pet=o$spet[o$sDates %in% o$rDates]
	o$ppet=o$spet[o$sDates %in% o$rPreDates]
	o$prain=o$srain[o$sDates %in% o$rPreDates]
	
	#calculate atributes
	attr(o$vi,'path')=dirname(o$vi[1])
	attr(o$rain,'path')=dirname(o$rain[1])
	attr(o$pet,'path')=dirname(o$pet[1])
	attr(o$ppet,'path')=dirname(o$ppet[1])
	attr(o$prain,'path')=dirname(o$prain[1])
	attr(o$vi,'files')=basename(o$vi)
	attr(o$rain,'files')=basename(o$rain)
	attr(o$pet,'files')=basename(o$pet)
	attr(o$ppet,'files')=basename(o$ppet)
	attr(o$prain,'files')=basename(o$prain)	
		
	#show
	showInfo(o)
	#return	
	o	
}


summarize=function(o){
	#parameters name in parent frame (to do a 'by reference')
	originalo=deparse(substitute(o))
	o$summarize=FALSE
	#outNames
	outNames=c('svi.stk','srain.stk','vimax','viMaxWhen')
	outNames=paste(o$pOut,'/',outNames,'.',o$driver,sep='')
	#reasterStack
	rasterStack(o$svi,outNames[1],interleave='BIP')
	rasterStack(o$srain,outNames[2],interleave='BIP')
	o$STKsvi=outNames[1]
	o$STKsrain=outNames[2]
	#viMax & viMaxWhen
	rgf.summary(o$vi,outNames[3],fun='MAX')
	aux=rgf.when(o$vi,outNames[3])
	writeGDAL(aux,outNames[4])
	o$viMax=outNames[3]
	o$viMaxWhen=outNames[4]
	#summary vi and rain series
	o$RESvi=rgf.summarize(o$vi)
	o$RESrain=rgf.summarize(o$rain)	
	o$summarize=TRUE
	assign(originalo,o,envir=parent.frame())
}

###############################################
# NAME: assessment
# PURPOSE:
# INPUTS:
# OUTPUTS:
assessment = function(o) {
	#parameters name in parent frame (to do a 'by reference')
	originalo=deparse(substitute(o))
	#set assessment flag to FALSE	
	o$rueEx=o$rueMed=o$aiEx=o$aiMed=''		
	o$assessment=FALSE
	assign(originalo,o,envir=parent.frame())
	#def files
	outNames=c('rueObsMe','rueObsEx','aiObsMe','aiObsEx')
	outNames=paste(o$pOut,'/',outNames,'.',o$driver,sep='')
	#try to make Indices	
	er=try({		
			rueMe=rueObsMe(o$rain,o$vi)
			writeGDAL(rueMe,outNames[1],drivername=o$driver,mvFlag=o$flag)
			iaMe=aiObsMe(o$rain,o$pet)
			writeGDAL(iaMe,outNames[3],drivername=o$driver,mvFlag=o$flag)
			rueEx=rueObsEx(o$rain,o$vi,o$prain,nMonths=o$acum)
			writeGDAL(rueEx,outNames[2],drivername=o$driver,mvFlag=o$flag)				
			iaEx=aiObsEx(o$rain,o$vi,o$pet,o$prain,o$ppet,nMonths=o$acum)	
			writeGDAL(iaEx,outNames[4],drivername=o$driver,mvFlag=o$flag)
	})
	#update o
	if (!class(er)=='try-error'){
		o$rueEx=outNames[1];
		o$rueMed=outNames[2];
		o$aiEx=outNames[3];
		o$aiMed=outNames[4];	
		o$assessment=TRUE 
		attr(o$assessment,'date')=format(Sys.time(),'%d %b %Y %H:%M:%S')
		aux=matrix(0, nrow = 4, ncol=7, dimnames = list(c('rueObsMe','rueObsEx','aiObsMe','aiObsEx'),c("Min","1st Qu","Median","Mean","3rd Qu","Max","NA's")))
		for (i in 1:4) {	
			aux1=summary(readGDAL(outNames[i],silent=TRUE)$band1)
			if (length(aux1)==6) aux1=c(aux1,0)
			aux[i,]=aux1
		}
		attr(o$assessment,'summary')=aux
		assign(originalo,o,envir=parent.frame())
	}
}



###############################################
# NAME: monitoring
# PURPOSE:
# INPUTS:
# OUTPUTS:
monitoring = function(o) {	
	#parameters name in parent frame (to do a 'by reference')
	originalo=deparse(substitute(o))
	#set monitoring flag to FALSE	
	o$f1=o$f2=o$f3=o$f4=o$f5=o$f6=o$f7=o$f8=o$f9=''			
	o$monitoring=FALSE
	assign(originalo,o,envir=parent.frame())
	#def files
	outNames=c('index','effect_time','effect_arid','veg_response','ta_single','tv_single','av_single')
	outNames=paste(o$pOut,'/',outNames,'.',o$driver,sep='')
	annualVis=paste(o$pOut,'/viMed',o$yIni:o$yEnd,'.',o$driver,sep='')
	annualIaMed=paste(o$pOut,'/aiMed',o$yIni:o$yEnd,'.',o$driver,sep='')
	annualTimes=paste(o$pOut,'/time',o$yIni:o$yEnd,'.',o$driver,sep='')
	#try to make Indices	
	er=try({	
		print('---------- Make annuals Vegetation Index Means')		
		rgf.summary(o$vi,annualVis,step=12,fun='MEAN',drivername=o$driver,mvFlag=o$flag)		
		#make aiObsMed by hidrologic years	
		print('---------- Make annuals Aridity Index')				
		n=o$rYears
		for (i in 0:(n-1)) {
			etp12=o$pet[(1:12)+i*12]
			rain12=o$rain[(1:12)+i*12]
			aiMe=aiObsMe(rain12,etp12,silent=T)
			writeGDAL(aiMe,annualIaMed[i+1],drivername=o$driver,mvFlag=o$flag)
		}		
		print('---------- Make annuals time series')
		#make annual time files		
		aux=readGDAL(annualVis[1],silent=T)	
		for (i in 0:(n-1)) {
			aux$band1=o$yIni+i
			writeGDAL(aux,annualTimes[i+1],drivername=o$driver,mvFlag=o$flag)
		}
		#make step by step regresion
		print('---------- Make Step by Step regresion')
		regStepRaster(annualVis,annualTimes,annualIaMed	,outNames,drivername=o$driver,mvFlag=o$flag)
	})
	#update o
	if (!class(er)=='try-error'){
		o$f1=outNames[1]
		o$f2=outNames[2]
		o$f3=outNames[3]
		o$f4=outNames[4]	
		o$f5=outNames[5]
		o$f6=outNames[6]
		o$f7=outNames[7]		
		o$monitoring=TRUE
		attr(o$monitoring,'date')=format(Sys.time(),'%d %b %Y %H:%M:%S')
		aux=matrix(0, nrow = 7, ncol=7, dimnames = list(c('index','effect_time','effect_arid','veg_response','ta_single','tv_single','av_single'),c("Min","1st Qu","Median","Mean","3rd Qu","Max","NA's")))
		for (i in 1:7) {
			aux1=summary(readGDAL(outNames[i],silent=TRUE)$band1)
			if (length(aux1)==6) aux1=c(aux1,0)
			aux[i,]=aux1
		}
		attr(o$monitoring,'summary')=aux		
		assign(originalo,o,envir=parent.frame())
	}
}


###############################################
# NAME: rueObsMe
# PURPOSE:
#     Reads n ndvi files and n rainfall grid files
#     Then Calculate the RueObsMe grid.
# INPUTS:
#       rainFl: File names list of rainfall grid 
#       viFl: File names list of vegetation index grid 
#       silent: logical Flag; if TRUE, comments outputs are supressed
# OUTPUTS:
#       Return RueObsMe as a SpatialGridDataFrame/SpatialPixelDataframe class. 
rueObsMe = function(rainFl, viFl, silent=FALSE) {
	
	if (length(rainFl)!=length(viFl)) stop('rainFl & viFl must have the same length')
	if (length(rainFl)%%12!=0) stop('rainFl must be multiply of 12')	
	if (length(rainFl)<12) stop('rainFl length must be greater or equal than 12')
	
	#Extraer el numero de elementos del array 
	n=length(rainFl) #num de meses en la serie
	nah=floor(n/12)    #num de years hidrologicos
	
	if (!silent) {
		print(paste('Procesing rueObsMe for', nah, 'years ---'))
		pb =txtProgressBar(min=0,max=n,char='*',width=20,style=3)
	}
	
	RueMed=0 
	
	# Para cada year hidrologico
	for (i in 0:(nah-1)) {	  
		SumRain=0
		SumNdvi=0
		#Para cada mes del year hidrologico
		for (j in 1:12){
			SumRain = SumRain+readGDAL(rainFl[i*12+j],silent=TRUE)$band1
			SumNdvi = SumNdvi+readGDAL(viFl[i*12+j],silent=TRUE)$band1
			if (!silent) setTxtProgressBar(pb, i*12+j)
		}
		SumNdvi=SumNdvi/12
		if (sum(SumRain == 0, na.rm=TRUE) > 0) SumRain=SumRain+(SumRain == 0) # cambio los pixels a 0mm de lluvia acumulada por 1mm
		RueMed=RueMed+(SumNdvi/SumRain)
	}
	RueMed=RueMed/nah
	aux=readGDAL(rainFl[1],silent=TRUE)
	aux$band1=RueMed
	if (!silent) close(pb)
	aux
}


##########################################################################
#
# NAME: GIaMed
#
# PURPOSE:
#     Reads n ndvi files and n rainfall files
#     Then Calculate the RUEMAX(ndvi_1..ndvi_n, rainfall_1..rainfall_n, m, errorfile) functiones with m#index
#
# INPUTS:
#       FilesArr: array [0..n#1] of String with name of IDRISI RST input files
#
#                 FilesArr[0]     # The no values mask filename
#                 FilesArr[1,12]  # The 12 precipitationes filename
#                 FilesArr[13,24] # The 12 evapotransipirationes files
#
# OUTPUTS:
#       Return Rue Maximun map,

aiObsMe=function (rainFl, petFl, FAO=FALSE, silent=FALSE){
	
	if (length(rainFl)!=length(petFl)) stop('rainFl & petFl must have the same length')
	if (length(rainFl)%%12!=0) stop('rainFl must be multiply of 12')	
	if (length(rainFl)<12) stop('rainFl length must be greater or equal than 12')
	
	#Extraer el numero de elementos del array 
	n=length(rainFl) #num de meses en la serie
	nah=floor(n/12) #num de years hidrologicos
	
	if (!silent){
		print(paste('Processing aiObsMe for', nah, 'years ---'))
		pb =txtProgressBar(min=0,max=n,char='*',width=20,style=3)
	}
	IaMed=0
	
	# Para cada year hidrologico
	for (i in 0:(nah-1)){
		SumRain=0
		SumPet=0
		#Para cada mes del year hidrologico
		for (j in 1:12){
			SumRain = SumRain+readGDAL(rainFl[i*12+j],silent=TRUE)$band1
			SumPet  = SumPet +readGDAL(petFl[i*12+j],silent=TRUE)$band1
			if (!silent) setTxtProgressBar(pb, i*12+j)
		}
		
		if (FAO) {
			if (sum(SumPet  == 0, na.rm=TRUE) > 0) SumPet =SumPet +(SumPet  == 0) # cambio los pixels a 0mm de Pet acumulada por 1mm
			IaMed=IaMed+(SumRain/SumPet)
		} else {
			if (sum(SumRain == 0, na.rm=TRUE) > 0) SumRain=SumRain+(SumRain == 0) # cambio los pixels a 0mm de lluvia acumulada por 1mm
			IaMed=IaMed+(SumPet/SumRain)
		}
	}
	
	IaMed=IaMed/nah
	aux=readGDAL(rainFl[1],silent=TRUE)
	aux$band1=IaMed
	if (!silent) close(pb)
	aux
}


############################333
# NAME: GRueMax
#
# PURPOSE:
#     Reads m ndvi files, m rainfall files, and a NMA map (numero de meses de acumulacion)
#     Then Calculate the RUEMAX
#
# INPUTS:
#       rainFl: File names list of rainfall grids 
#       viFl: File names list of vegetation index grids
#       preRainFl:  File names list of previous rainfall grids 
#       Nma: cte de acumulacion
#       silent: logical Flag; if TRUE, comments outputs are supressed
# OUTPUTS:
#       Return Rue Maximun map, and the Rue Max Month map
# TODO: Fallo con acum=1 .. ver = en aiObsEx

rueObsEx = function (rainFl, viFl, preRainFl, nMonths=6, silent=FALSE){
	
	if (length(rainFl)!=length(viFl)) stop('rainFl & viFl must have the same length')
	if (length(rainFl)%%12!=0) stop('rainFl must be multiply of 12')	
	if (length(rainFl)<12) stop('rainFl length must be greater or equal than 12')
	if (is.character(nMonths)) {
		faux=readGDAL(nMonths,silent=TRUE)
		m=faux$band1
		nMonths=max(m)
	} else {
		faux=readGDAL(rainFl[1],silent=TRUE)
		faux$band1=nMonths
		m=faux$band1
	}
	
	if (length(preRainFl)<nMonths) stop(paste('preRain length mus be at least ',nMonths,'elements'))
	
	#Extraer el numero de elementos del array 
	n=length(rainFl) #num de meses en la serie
	nah=floor(n/12) #num de years hidrologicos
	
	RainFifo=0 
	NdviMax=0
	NdviMaxMon=0
	Aux=0
	RueMax=0    
	
	# Calcula NdviMax y NdviMaxMonth
	if (!silent) {
		print(paste('Processing rueObsEx for', nah, 'years ---'))
		print('Processing vegetation index files')
		pb =txtProgressBar(min=1,max=n,char='*',width=20,style=3)
	} 
	for (i in 1:n) {
		Aux=readGDAL(viFl[i],silent=TRUE)$band1
		NdviMax=pmax(Aux,NdviMax)
		Msk=(Aux == NdviMax)
		NdviMaxMon=NdviMaxMon*(!Msk)+i*Msk #Mes indicado de 0 a n-1
		if (!silent) setTxtProgressBar(pb, i)
	}
	 
	if (!silent) {
		close(pb)
		print('Processing rain files')
		pb =txtProgressBar(min=1,max=nMonths,char='*',width=20,style=3)
	}
	#relleno fifo con los nMonths primeros meses    
	for (i in 1:nMonths) {
		RainFifo=cbind(RainFifo,readGDAL(preRainFl[i],silent=TRUE)$band1)
		setTxtProgressBar(pb, i)
	}
	RainFifo=RainFifo[,-1] #eliminamos primera columna de 0s
	
	if (!silent) {
		close(pb)			
		pb =txtProgressBar(min=1,max=n,char='*',width=20,style=3)
	}
#browser()
	#for each month in the serie
	for (i in 1:n) {
		#Calculate the mask of NdviMaxMon for month i
		#pixeles donde el ndvi maximo se obtuvo en el mes numero i
		MskMonth=(NdviMaxMon == i)       
		# Para cada valor de los nMonths posibles
		SumRain=0
		for (j in 1:nMonths) {
			#Calcular imagenes de lluvia acumulada en j meses
			SumRain=SumRain+RainFifo[,nMonths+1-j] #lo acumulamos
			#Calcula mascara de pixeles con m=j
			MsknMonths=(m == j)#
			if (sum(MsknMonths)>0) {
				#Calcular RueMax
				if (sum(SumRain == 0, na.rm=TRUE) > 0) SumRain=SumRain+(SumRain == 0) # cambio los pixels a 0mm de lluvia acumulada por 1mm
				RueMax=RueMax*(!MskMonth)+(NdviMax/SumRain)*MsknMonths*MskMonth #overlay RueMax
			}
		}
		#Desplazar Fifo e insertar mes
		RainFifo=cbind(RainFifo[,-1],readGDAL(rainFl[i],silent=TRUE)$band1)
		if (!silent) setTxtProgressBar(pb, i)
	}
	
	faux$band1=RueMax
	if (!silent) close(pb)
	faux
}

# NAME:
#       GIAMax
#
# PURPOSE:
#     Reads n ndvi files and n rainfall files
#     Then Calculate the RUEMAX(ndvi_1..ndvi_n, rainfall_1..rainfall_n, m, errorfile) functiones with m-index
#
# INPUTS:
#       FilesArr: array [0..n-1] of String with name of IDRISI RST input files
#
#                 FilesArr[0]     - The no values mask filename
#                 FilesArr[1,12]  - The 12 precipitationes filename
#                 FilesArr[13,24] - The 12 evapotransipirationes files
#
# OUTPUTS:
#       Return Rue Maximun map,

aiObsEx = function (rainFl, viFl, petFl, preRainFl, prePetFl, FAO=FALSE, nMonths=6, silent=FALSE) {
	if (length(rainFl)!=length(viFl)) stop('rainFl & viFl must have the same length')
	if (length(rainFl)%%12!=0) stop('rainFl must be multiply of 12')	
	if (length(rainFl)<12) stop('rainFl length must be greater or equal than 12')
	if (is.character(nMonths)) {
		faux=readGDAL(nMonths,silent=TRUE)
		m=faux$band1
		nMonths=max(m)
	} else {
		faux=readGDAL(rainFl[1],silent=TRUE)
		faux$band1=nMonths
		m=faux$band1
	}
	
	if (length(preRainFl)<nMonths) stop(paste('preRain length mus be at least ',nMonths,'elements'))
	
	#Extraer el numero de elementos del array 
	n=length(rainFl) #num de meses en la serie
	nah=floor(n/12) #num de years hidrologicos
	
	RainFifo=0 
	PetFifo=0  
	viMax=0
	viMaxMon=0
	Aux=0
	IaMax=0
	
	# Calcula NdviMax y NdviMaxMonth
	if (!silent) {
		print(paste('Processing aiObsEx for', nah, 'years ---'))
		print('Processing vegetation index files')
		pb =txtProgressBar(min=1,max=n,char='*',width=20,style=3)
	}
	for (i in 1:n) {
		Aux=readGDAL(viFl[i],silent=TRUE)$band1
		viMax=pmax(Aux,viMax)
		Msk=(Aux == viMax)
		viMaxMon=viMaxMon*(!Msk)+i*Msk #Mes indicado de 0 a n-1
		if (!silent) setTxtProgressBar(pb, i)
	}
	if (!silent) {
		close(pb)
		print('Processing data')
	}
	pb =txtProgressBar(min=1,max=nMonths,char='*',width=20,style=3)
	#relleno fifo con los nMonths primeros meses    
	for (i in 1:nMonths) {
		RainFifo=cbind(RainFifo,readGDAL(preRainFl[i],silent=TRUE)$band1)
		PetFifo=cbind(PetFifo,readGDAL(prePetFl[i],silent=TRUE)$band1)
		if (!silent) setTxtProgressBar(pb, i)
	}
	RainFifo=RainFifo[,-1] #eliminamos primera columna de 0s
	PetFifo=PetFifo[,-1] #eliminamos primera columna de 0s
	
	if (!silent) {
		close(pb)		
		pb =txtProgressBar(min=1,max=n,char='*',width=20,style=3)
	}
	#for each month in the serie
	for (i in 1:n) {
		#Calculate the mask of viMaxMon for month i
		#pixeles donde el vi maximo se obtuvo en el mes numero i
		MskMonth=(viMaxMon == i)       
		# Para cada valor de nMonths de los 12 posibles
		SumRain=0
		SumPet=0
		for (j in 1:nMonths) {
			#Calcular imagenes de lluvia acumulada en j meses
			SumRain=SumRain+RainFifo[,nMonths+1-j] #lo acumulamos
			SumPet=SumPet+PetFifo[,nMonths+1-j] #lo acumulamos
			#Calcula mascara de pixeles con nMonths=j
			MsknMonths=(m == j)#
			if (sum(MsknMonths)>0) {
				#Calcular IaMax
				if (FAO) {
					if (sum(SumPet  == 0, na.rm=TRUE) > 0) SumPet =SumPet +(SumPet  == 0) # cambio los pixels a 0mm de Pet acumulada por 1mm
					IaMax=IaMax*(!MskMonth)+(SumRain/SumPet)*MsknMonths*MskMonth #overlay RueMax
				} else {
					if (sum(SumRain == 0, na.rm=TRUE) > 0) SumRain=SumRain+(SumRain == 0) # cambio los pixels a 0mm de lluvia acumulada por 1mm
					IaMax=IaMax*(!MskMonth)+(SumPet/SumRain)*MsknMonths*MskMonth #overlay RueMax
				}
			}
		}
		#Desplazar Fifo e insertar mes
		RainFifo=cbind(RainFifo[,-1],readGDAL(rainFl[i],silent=TRUE)$band1)
		PetFifo=cbind(PetFifo[,-1],readGDAL(petFl[i],silent=TRUE)$band1)
		if (!silent) setTxtProgressBar(pb, i)
	}
	faux$band1=IaMax
	if (!silent) close(pb)
	faux
}

###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################
rasterStack=function(inFl,outFN,asc=FALSE,zip=FALSE,dec=3,interleave='BIP',silent=FALSE){
	
	#comprueba condiciones de error
	if (length(inFl)<2) stop('inFl must be at least of length 2')
	if (nchar(outFN)==0) stop('Empty output file name')
	if (interleave %in% c('BIL','BIP','BSQ') == FALSE) stop('Invalid interleave. Valid are BIL, BIP, BSQ')
	if ((interleave %in% c('BIL','BSQ')) & (!missing(asc) || !missing(zip) || !missing(dec))) { 
		asc=FALSE
		zip=FALSE
		dec=3
		print('asc,zip,dec options are ignored in interleave mode BIL or BSQ')
	}
	
	#num de imagenes en la lista 
	nimg=length(inFl)
	
	#num de filas y columnas de la imagen
	rows=GDALinfo(inFl[1])[1]
	cols=GDALinfo(inFl[1])[2]
	
	#calculo size del buffer de lectura para que lea bloques de 100MB
	linesToRead=ceiling(100000000/(cols*nimg*8))
	
	#num de bloques de size linesToRead en la imagen
	nblocks=ceiling(rows/linesToRead)
	
	#inicio progressbar
	if (!silent) pb=txtProgressBar(min=0,max=nblocks*nimg,char='*',width=20,style=3)
	
	#abrir fichero de salida en modo write + (texto, texto comprimido o binario)
	if (asc) { 
		if (zip) {
			ft=gzfile(outFN,'wt')
		} else ft=file(outFN,'wt')
	} else ft=file(outFN,'wb')
	
	if (interleave=='BSQ') {
		stop('sorry, not implemented yet...')
	}  else {
		
		#leemos lineToRead (ej. 200 lineas) de cada imagen 
		#anexamos cada bloque leido a la columna de un dataframe 
		#escribimos en disco el datafame con la opcion de anexar 
		#el proceso se repite para los nblocks necesarios
		for (j in 0:(nblocks-1)){
			#iniciamos dataframe de salida
			outdf=0
			for (i in 1:nimg){
				aux=GDAL.open(inFl[i],TRUE)  
				#offsets
				offsetIni=j*linesToRead+1
				offsetFin=(j+1)*linesToRead
				#no rebasar fin de archivo
				if (offsetFin>rows) offsetFin=rows
				#lee bloque de la imagen SIN cargarla completamente en memoria 
				Strip=aux[offsetIni:offsetFin,]
				#cerrar conexion
				GDAL.close(aux)
				#anexar como columna al dataframe de salida    
				outdf=cbind(outdf,Strip$band1)  
				#actualizo progressbar
				if (!silent) setTxtProgressBar(pb, i+nimg*j)
			} #for
			#convertir outdf a vector, y darle la forma adecuada para guardarlo en disco
			outdf=outdf[,-1]
			if (asc) {
				outdf=round(outdf,dec)
				write.table(outdf,ft,row.names=FALSE,col.names=FALSE,sep=',')
			} else {
				if (interleave=='BIL') writeBin(as.vector(as.matrix(outdf)),ft,size=4)
				if (interleave=='BIP' & !asc) writeBin(as.vector(t(as.matrix(outdf))),ft,size=4)
			}
		} #for
		close(ft)
		close(pb)
	}
}

###############################################
# NAME: 
# PURPOSE:
# INPUTS:regStepRaster(vi,y,ia,fl,drivername='RST',mvFlag=-1)
# OUTPUTS:
###############################################
regStepRaster=function(ndviFl,timeFl,aridFl,outFl,silent=FALSE,...){
	#comprueba condiciones de error
	if (length(outFl)!=7) stop('outFl may be a list of seven filenames')
	if (2*length(ndviFl)!=(length(timeFl)+length(aridFl))) stop('ndviFl, tempFl and aridFl, should have equal length')
	
	#Stack by pixel the filelists
	tmpFn1=tempfile()
	tmpFn2=tempfile()
	tmpFn3=tempfile()
	tmpFn4=tempfile()	#fichero temporal para almacenar resultados de la regresion paso a paso
	rasterStack(ndviFl,tmpFn1,interleave='BIP')
	rasterStack(timeFl,tmpFn2,interleave='BIP')
	rasterStack(aridFl,tmpFn3,interleave='BIP')
	
	#image info
	bands=length(ndviFl)	
	rows=GDALinfo(ndviFl[1])[1]
	cols=GDALinfo(ndviFl[1])[2]
	#browser()
	#calculo size del buffer de lectura para que lea bloques de 5000 elementos aproximadamente
	items=bands*rows*cols
	#este es un size optimo para la funcion 'by'
	pixelsToRead=trunc(5000/bands)	
	#num de bloques de size itemsToRead
	nblocks=trunc(items/(pixelsToRead*bands))
	#tamaño del ultimo bloque
	rest=rows*cols-nblocks*pixelsToRead
	
	if (!silent) pb=txtProgressBar(min=0,max=nblocks,char='*',width=20,style=3)
	
	depf=file(tmpFn1,'rb')
	in1f=file(tmpFn2,'rb')	
	in2f=file(tmpFn3,'rb')
	outf=file(tmpFn4,'w')
	
#	cat(sprintf('\n%dx%dx%d, pixels=%d,   items=%d',rows,cols,bands,rows*cols,rows*cols*bands))#
	#cat(sprintf('\n 	  pixToRead=%d, bloques=%d, resto=% dpixels',pixelsToRead,nblocks,rest))
	#cat('\n')
	
	#por cada (bloque + 1) 
	for (i in 0:nblocks) {
		if (i == nblocks) {pixelsToRead=rest}
		#leer un linestoread de lineas del fichero de entrada
		Y=readBin(depf,numeric(),pixelsToRead*bands,size=4)
		X1=readBin(in1f,numeric(),pixelsToRead*bands,size=4)
		X2=readBin(in2f,numeric(),pixelsToRead*bands,size=4)
	
		df=cbind(Y,X1,X2,pixel=rep(1:pixelsToRead,each=bands))
		rm(Y,X1,X2)
	
		#calcular regresion multiple
		cn=regStepDF(df)
		
		#escribir salida
		write.table(round(cn,4),append=TRUE,sep='\t',file=outf,col.names=FALSE,row.names=FALSE)
	
		#actualizo progressbar
		if (!silent) setTxtProgressBar(pb,i)
	}
	
	close(depf)
	close(in1f)
	close(in2f)
	flush(outf)
	close(outf)
	
	if (!silent) close(pb)

	aux1=read.table(tmpFn4,header=FALSE,sep='\t')
	aux2=readGDAL(ndviFl[1],silent=TRUE)
	print('writing output files')
	for (i in 1:7) {
		aux2$band1=aux1[,i]
		writeGDAL(aux2,outFl[i],...)
	}	
	file.remove(tmpFn1)	
	file.remove(tmpFn2)	
	file.remove(tmpFn3)	
	file.remove(tmpFn4)	
}

###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################
regStepDF=function (X){
	dimX=dim(X)[1]
	cols=max(X[,4])
	ocases=length(unique(X[,4]))
	#quitar casos con algun NA en alguna de las bandas
	aux=unique(X[is.na(X[,1]*X[,2]*X[,3]),4]) #lista de pixel con algun dato a NA
	X[X[,4] %in% aux,1]=NA	
	X=na.omit(X) 

	ndimX=dim(X)[1]
	index=unique(X[,4])
	#si todo el dataframe es iniutilizable
	if (ndimX==0) return(matrix(NA,ocases,7))
	
	#num de bandas en la matriz
	N=sum(X[,4]==X[,4][1])
	ncases=length(index)
	
	
	RX=0
	#print(paste('dimX:',dimX,' dimXna:',ndimX,' ncases:',ncases))
	#browser()
	# Correlaciones simples entre variables segun Box 15.2
	#aux=(unlist(by(X[,1:3],X[,4],cor)))
	for (i in 0:(ncases-1)){		
		RX=rbind(RX,cor(X[(i*N+1):((i+1)*N),1:3])[c(6,2,3)])
	}
	RX=RX[-1,]
	
	#caso especial...
	#si en X hay un solo pixel valido, y el resto son NAs, 
	#RX adopta la forma [ 0, 0, 0] 
	#                   [ A, B, C]
	#y al quitar la primera linea pasa de ser una matriz a ser un vector...
	#así que hay que forzar a ser matriz
	if (class(RX)=='numeric') {RX= as.matrix(t(RX))}	
	
	#dim(aux)=c(9,ncases)
	#RX=t(aux[c(6,2,3),])
	#colnames(RX)=c('Rx1x2','Rx1Y','Rx2Y')
	
	# Coeficientes standard de regresion parcial
	# con ecuacion YP=BPY1*X1P+BPY2*X2P
	
	MBPY1=(RX[,2]-RX[,3]*RX[,1])/(1-RX[,1]^2)
	MBPY2=(RX[,3]-RX[,2]*RX[,1])/(1-RX[,1]^2)
	
	# Estadisticos de significacion de correlacion simple segun Box 15.4 con N<50
	DFS = N - 2
	MTS = RX * sqrt(DFS / (1 - RX^2))
	
	# Coeficiente de determinacion multiple
	MR2Y12 = RX[,2] * MBPY1 + RX[,3] * MBPY2
	
	# Estadistico de significacion de la regresion multiple con 2 variables independientes
	DFNUM = 2
	DFNUM2= 1
	DFDEN = N - 3
	MF = (MR2Y12 / 2) / ((1 - MR2Y12) / DFDEN)
	
	# Incremento en determinacion al meter la segunda variable
	ARX=abs(RX)
	MRXY=ifelse(ARX[,2]>ARX[,3],RX[,2],RX[,3])
	MV2 =ifelse(ARX[,2]>ARX[,3],2,1)
	MF2 = (MR2Y12 - MRXY^2) / ((1 - MR2Y12) / DFDEN)
	
	# Distribuciones de Probabilidad
	MPrR2Y12 = 1-pf(MF,DFNUM,DFDEN)
	MPrR2Y12i= 1-pf(MF2,DFNUM2,DFDEN)
	
	MPr = 2*(1-pt(abs(MTS),DFS))
	
	#MPrX1X2 = 2*(1-pt(abs(MTS[,1]),DFS))
	#MPrX1Y  = 2*(1-pt(abs(MTS[,2]),DFS))
	#MPrX2Y  = 2*(1-pt(abs(MTS[,3]),DFS))
	
	#matriz de salida del efecto de la aridez y el tiempo q[,1] tiempo, q[,2] aridez, q[,3] respuesta de la vegetacion (1,2,3 o 4)
	q=matrix(-1,ncases,7) #casos sin NA
	nq=matrix(NA,cols-ncases,7) #casos con NA
	colnames(q) =c('index','effect_time','effect_arid','veg_response','ta_single','tv_single','av_single')
	colnames(nq)=c('index','effect_time','effect_arid','veg_response','ta_single','tv_single','av_single')
	
	#indetifica el pixel al que peretenece el resultado
	q[,1]=index
	nq[,1]=(1:cols)[-index]
	# wi almacena los indices de los casos donde se cumple la condicion...
	# 1a condicion 
	wi=which((MPrR2Y12<=0.1) & (MPrR2Y12i<=0.1))
	q[wi,2:3]=cbind(MBPY1[wi],MBPY2[wi])
	# 2a condicion 
	wi=c(which((MPrR2Y12<=0.1) & (MPrR2Y12i>0.1) & (MV2==2)),which((MPrR2Y12>0.1) & (MPr[,2]<=0.1)))
	q[wi,2:3]=cbind(RX[wi,2],0)
	# 3a condicion 
	wi=c(which((MPrR2Y12<=0.1) & (MPrR2Y12i>0.1) & (MV2==1)),which((MPrR2Y12>0.1) & (MPr[,2]>0.1) & (MPr[,3]<=0.1)))
	q[wi,2:3]=cbind(0,RX[wi,3])
	# 4a condicion 
	wi=which(((MPrR2Y12>0.1) & (MPr[,2]>0.1) & (MPr[,3]>0.1)))
	q[wi,2:3]=cbind(0,0)
	# calculo de la respuesta de la vegetacion (rv) en funcion del efecto del tiempo (et) y la aridez (ea)
	#  et  ea  rv
	# no 0  0-> 1
	# no 0 no 0 -> 3
	#  0 no 0 -> 2
	#  0  0-> 4
	
	q[,4]=ifelse(q[,2]!=0,ifelse(q[,3]==0,1,3),ifelse(q[,3]==0,4,2))
	
	# calculo de variables mostrando el efecto simple
	q[,5]=ifelse(MPr[,1]<=0.1,RX[,1],0)
	q[,6]=ifelse(MPr[,2]<=0.1,RX[,2],0)
	q[,7]=ifelse(MPr[,3]<=0.1,RX[,3],0)
	
	q=rbind(q,nq)
	q=q[order(q[,1]),]
	# ----------vector de resultados
	return(q)
	
	# ----------devuelve vector completo para debug
	#round(c(R12, R1Y, R2Y, DFS, TS12, TS1Y, TS2Y, PrX1X2,PrX1Y,PrX2Y,BPY1, BPY2, BY1, BY2, A, R2Y12, DFNUM, DFDEN, F,PrR2Y12, V2, DFNUM2, F2,PrR2Y12i, DFM, TM1, TM2, TMP1, TMP2),6)
}
