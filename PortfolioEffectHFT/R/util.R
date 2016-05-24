util_getComputeTime<-function(time=c("timeMax","timeLeft")){
	util_validate()
	clientConnection=getOption('clientConnection')
	result=.jcall(clientConnection,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getComputeTimeLeft")
	util_checkErrors(result)
	result=.jcall(result,returnSig="S", method="getInfoParam",time[1])
	return(result)
}

util_POSIXTimeToDate<-function(time){
	if(is.numeric(time)){
	DateTimeUtil=.jnew("com.portfolioeffect.quant.client.util.DateTimeUtil")
	time=.jcall(DateTimeUtil,returnSig='[S',method="POSIXTimeToDateStr",.jlong(time))
		}
	time
}

util_dateToPOSIXTime<-function(time){
	if(!is.numeric(time)){
	DateTimeUtil=.jnew("com.portfolioeffect.quant.client.util.DateTimeUtil")
	check=1414818000000-.jcall(DateTimeUtil,returnSig='[J',method="toPOSIXTime",'2014-11-01 01:00:00')
	time=.jcall(DateTimeUtil,returnSig='[J',method="toPOSIXTime",time)+check
}
time
}

util_screenshot<-function(path, delaySec=15){
	Sys.sleep(delaySec)
	J("com.portfolioeffect.quant.client.util.Screenshots")$takeScreenshot(path)
}

util_setCredentials<-function(username,password,apiKey,hostname="snowfall04.snowfallsystems.com"){
	way<-switch(Sys.info()[['sysname']],
			Windows= {paste(Sys.getenv("APPDATA"),"\\ice9",sep="")},
			Linux  = {paste(Sys.getenv("HOME"),"/.ice9",sep="")},
			Darwin = {paste(Sys.getenv("HOME"),"/.ice9",sep="")})
	way<-gsub("\\","/",way,fixed = T)
	way<-paste0(way, "/data")
	if (!file.exists(way)){
		dir.create(way,recursive=TRUE)
	}
	unlink(paste0(way,"/login.RData"), recursive = TRUE, force = FALSE)
	login<-c(username,password,apiKey,hostname)
	save(login, file=paste0(way, "/login.RData"))
	rm(login)
	if(!is.null(options('clientConnection')$clientConnection)){
		clientConnection=getOption('clientConnection')
		if(clientConnection == NULL){
			options("clientConnection"=.jnew("com.portfolioeffect.quant.client.ClientConnection"))
			is.null(options('clientConnection')$clientConnection)
		}
		.jcall(clientConnection,returnSig="V", method="setUsername",username)
		.jcall(clientConnection,returnSig="V", method="setPassword",password)
		.jcall(clientConnection,returnSig="V", method="setApiKey",apiKey)
		.jcall(clientConnection,returnSig="V", method="setHost",hostname)
		resultTemp=.jcall(clientConnection,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="restart")
		util_checkErrors(resultTemp)
	}
	util_validate()
	clientConnection=getOption('clientConnection')
	resultTemp=.jcall(clientConnection,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="restart")
	util_checkErrors(resultTemp)
}

util_cleanCredentials<-function(){
	way<-switch(Sys.info()[['sysname']],
			Windows= {paste(Sys.getenv("APPDATA"),"\\ice9",sep="")},
			Linux  = {paste(Sys.getenv("HOME"),"/.ice9",sep="")},
			Darwin = {paste(Sys.getenv("HOME"),"/.ice9",sep="")})
	way<-gsub("\\","/",way,fixed = T)
	way<-paste0(way, "/data")
	if (!file.exists(way)){
		dir.create(way,recursive=TRUE)
	}	
	unlink(paste0(way,"/login.RData"), recursive = TRUE, force = FALSE)
}

util_validate<- function(argList=NULL) {
	names<-names(argList)
	for(name in names){
		if((name %in% c('portfolio'))&(!is(argList[[name]],'portfolio'))){
			stopMessage("NOT_PORTFOLIO_CLASS")
		}
		if((name %in% c('optimizer'))&(!is(argList[[name]],'optimizer'))){
			stopMessage("NOT_OPTIMIZER_CLASS")
		}
		if((name %in% c('symbol','symbol1','symbol2','metricType','constraintType','constraintMertic','goal','direction'))&(!is(argList[[name]],'character'))){
			stopMessage("OBJECT_NOT_CHARACTER_CLASS",names=name)
		}
		if((name %in% c('lag','confidenceInterval','thresholdReturn','confidenceIntervalA','confidenceIntervalB','pValueLeft','pValueRight','nPoints','portfolioValue'))&(!is(argList[[name]],'numeric'))){
			stopMessage("OBJECT_NOT_NUMERIC_CLASS",names=name)
		}
		if((name %in% c('lag','confidenceInterval','thresholdReturn','confidenceIntervalA','confidenceIntervalB','pValueLeft','pValueRight','nPoints','portfolioValue'))&(!(length(argList[[name]])==1))){
			stopMessage("OBJECT_NOT_SINGLE_NUMBER",names=name)
		}
		if(name %in% c('lag','confidenceInterval','thresholdReturn','confidenceIntervalA','confidenceIntervalB','pValueLeft','pValueRight','nPoints')){
			if((argList[[name]]<0)){
			stopMessage("OBJECT_NOT_POSITIVE_NUMBER",names=name)
		}
	}
	}
	way<-switch(Sys.info()[['sysname']],
			Windows= {paste(Sys.getenv("APPDATA"),"\\ice9",sep="")},
			Linux  = {paste(Sys.getenv("HOME"),"/.ice9",sep="")},
			Darwin = {paste(Sys.getenv("HOME"),"/.ice9",sep="")})
	way<-gsub("\\","/",way,fixed = T)
	way<-paste0(way, "/data/login.RData")
	if(!file.exists(way)){
		stopMessage('FILE_CREDENTIALS_NO_EXISTS')
	}
	
	if(is.null(options('clientConnection')$clientConnection)){
		options("clientConnection"=.jnew("com.portfolioeffect.quant.client.ClientConnection"))
		clientConnection=getOption('clientConnection')
		login<-NULL;
		load(way)
		.jcall(clientConnection,returnSig="V", method="setUsername",login[1])
		.jcall(clientConnection,returnSig="V", method="setPassword",login[2])
		.jcall(clientConnection,returnSig="V", method="setApiKey",login[3])
		.jcall(clientConnection,returnSig="V", method="setHost",login[4])
		
	}else{
		clientConnection=getOption('clientConnection')
		if(clientConnection==NULL){
		options("clientConnection"=.jnew("com.portfolioeffect.quant.client.ClientConnection"))
		clientConnection=getOption('clientConnection')
			login<-NULL;
			load(way)
			.jcall(clientConnection,returnSig="V", method="setUsername",login[1])
			.jcall(clientConnection,returnSig="V", method="setPassword",login[2])
			.jcall(clientConnection,returnSig="V", method="setApiKey",login[3])
			.jcall(clientConnection,returnSig="V", method="setHost",login[4])
		}
		else {
			login<-NULL;
			load(way)
			.jcall(clientConnection,returnSig="V", method="setUsername",login[1])
			.jcall(clientConnection,returnSig="V", method="setPassword",login[2])
			.jcall(clientConnection,returnSig="V", method="setApiKey",login[3])
			.jcall(clientConnection,returnSig="V", method="setHost",login[4])
		}
	}
}

util_checkErrors<-function(result){
	if(getErrorStatus(result)){
		Message=getErrorMessage(result)
		n=nchar(Message)
		c=paste(array("#",dim=min(abs(round(((n-15)/2)+0.01))),0),collapse="")
		cc=paste(array("#",dim=min(n+(n-15)%%2),0),collapse="")
		k=paste(c," ERROR MESSAGE ",c,sep="")
		stop(paste("",k,Message,cc,sep="\n"),call.=FALSE)
	}
}

printStatus<-function(){
	clientConnection=getOption('clientConnection')
	temp<-.jcall(clientConnection,returnSig="S", method="getStatus")
	if(temp!=""){print(temp)}
}

getResult<-function(data){
	util_checkErrors(data)
	dataNames<-.jcall(data,returnSig="[S", method="getDataNames")
	result=NULL
	for(dataName in dataNames){
		dataType<-.jcall(data,returnSig="S", method="getDataType",dataName)
		resultTemp<-switch(dataType,
#		NULL          =,
		DOUBLE =.jcall(data,returnSig="D", method="getDouble", dataName),
		DOUBLE_VECTOR =.jcall(data,returnSig="[D", method="getDoubleArray", dataName),
		DOUBLE_MATRIX =.jcall(data,returnSig="[[D", method="getDoubleMatrix", dataName, simplify=TRUE),
		INT_VECTOR    =.jcall(data,returnSig="[I", method="getIntArray", dataName),
#		INT_MATRIX    =.jcall(data,returnSig="[[I", method="getDoubleMatrix", dataName, simplify=TRUE),
		LONG_VECTOR   =  .jcall(data,returnSig="[J", method="getLongArray", dataName),
#		LONG_MATRIX   =.jcall(data,returnSig="[[J", method="getDoubleMatrix", dataName, simplify=TRUE),
#		FLOAT_VECTOR  =.jcall(data,returnSig="[D", method="getFloatArray", dataName),
#		FLOAT_MATRIX  =.jcall(data,returnSig="[[D", method="getDoubleMatrix", dataName, simplify=TRUE),
#		STRING        =.jcall(data,returnSig="[D", method="getDoubleArray", dataName),
		STRING_VECTOR =.jcall(data,returnSig="[S", method="getStringArray", dataName),
		PORTFOLIO     =.jcall(data,returnSig="Lcom/portfolioeffect/quant/client/portfolio/Portfolio;", method="getPortfolio", "portfolio"))
if(dataType=="PORTFOLIO"){
	result=resultTemp
}else{
result=cbind(result,resultTemp)	
}
}
if(NROW(result)>1){
colnames(result)<-dataNames
}
return(result)
}

getErrorStatus<-function(result){
	.jcall(result,returnSig="Z", method="hasError")
}
getErrorMessage<-function(result){
	.jcall(result,returnSig="S", method="getErrorMessage")
}
getTimeMilliSec<-function(result){
	.jcall(result,returnSig="[J", method="getLongArray", "time")
}
getResultValuesLong<-function(result){
	.jcall(result,returnSig="J", method="getLastLong", "value")
}
getResultValuesDoubleArray<-function(result){
.jcall(result,returnSig="[D", method="getDoubleArray", "value")
}
getResultValuesDoubleArrayWithTime<-function(result){
	result<-cbind(.jcall(result,returnSig="[J", method="getLongArray", "time"),.jcall(result,returnSig="[D", method="getDoubleArray", "value"))
	colnames(result)<-c("Time","Value")
	return(result)
}
getResultValuesDouble2DArray<-function(portfolio,result){
	result<-cbind(.jcall(result,returnSig="[J", method="getLongArray", "time"),.jcall(result,returnSig="[[D", method="getDoubleMatrix", "value", simplify=TRUE))
	colnames(result)<-c("Time",portfolio_symbols(portfolio))
	return(result)

}
getResultValuesPortfolio<-function(result){
	.jcall(result,returnSig="Lcom/portfolioeffect/quant/client/portfolio/Portfolio;", method="getPortfolio", "portfolio")
}

