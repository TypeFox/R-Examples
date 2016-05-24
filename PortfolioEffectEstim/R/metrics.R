setClass("estimator",
		slots = c(java="jobjRef")
)
estimator_mainSettings<-function(estimator){
	estimator_settings(estimator,
			jumpsModel = "moments",
			resultsSamplingInterval = "1s",
			inputSamplingInterval="none")
}
estimator_settings<-function(estimator,...){
	util_validate()
	isList=FALSE
	try({if(is.list(...)){isList=TRUE}},silent = TRUE)
	if(isList){change=(...)}else{change=list(...)}
	removeList<-NULL
	for(i in 1:length(change)){
		if(!(names(change[i]) %in% c( "jumpsModel","resultsSamplingInterval","inputSamplingInterval"))){
			stopMessage('WRONG_SETTINGS_ARGUMENTS')
		}
		if(names(change[i]) %in% c( 'resultsSamplingInterval','inputSamplingInterval')){
			switch(names(change[i]),
					resultsSamplingInterval = {
						if(is.character(change[[i]])){
							change$samplingInterval=change[[i]]
						}else{
							change$samplingInterval=.jlong(change[[i]])
						}},
					inputSamplingInterval  = {
						change$priceSamplingInterval=change[[i]]
					}
			)
			removeList<-c(removeList,i)
		}
	}
	change[removeList]<-NULL
	for(i in 1:length(change)){
		.jcall(estimator@java,returnSig="V", method="setParam",names(change[i]),paste(change[[i]]))
	}
}

estimator_getSettings<-function(estimator){
	
	temp<-list()
	temp$jumpsModel<-.jcall(estimator@java,returnSig="S", method="getParam","jumpsModel")
	temp$resultsSamplingInterval<-.jcall(estimator@java,returnSig="S", method="getParam","samplingInterval")
	temp$inputSamplingInterval<-.jcall(estimator@java,returnSig="S", method="getParam","priceSamplingInterval")
	temp
}
#
# Rolling integrated variance
#
estimator_create<-function(asset=NULL,fromTime=NULL,toTime=NULL,priceData=NULL){
	util_validate()
	clientConnection=getOption('clientConnection')
	estimator=new("estimator", java=.jnew("com.portfolioeffect.quant.client.portfolio.Estimator",clientConnection))
	
	
	if(is.null(asset) & is.null(fromTime)& is.null(toTime)& is.null(priceData)){
		stop('No arguments provided, please check required arguments list.')
	}
	if(((class(asset)=="matrix")&(is.null(priceData)))||((class(priceData)=="matrix")&(is.null(asset)))){
		if((class(asset)=="matrix")&(is.null(priceData))){
			priceData=asset
		}
		util_validate()
		clientConnection=getOption('clientConnection')
		estimator=new("estimator", java=.jnew("com.portfolioeffect.quant.client.portfolio.Estimator",clientConnection))
		result<-.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addAsset", as.double(priceData[,2]),.jlong(priceData[,1]))
		util_checkErrors(result)
		estimator_mainSettings(estimator)
		return(estimator)
	}
	if((class(fromTime)=="character")&(class(toTime)=="character")&(is.null(priceData))){
		
		util_validate()
		clientConnection=getOption('clientConnection')
		estimator=new("estimator", java=.jnew("com.portfolioeffect.quant.client.portfolio.Estimator",clientConnection))
		result<-.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setFromTime", fromTime)
		util_checkErrors(result)
		result<-.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setToTime", toTime)
		util_checkErrors(result)
		result<-.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addAsset", asset)
		util_checkErrors(result)
		estimator_mainSettings(estimator)
		return(estimator)
	}
	if((class(asset)=="estimator")&(is.null(fromTime))&(is.null(toTime))&(is.null(priceData))){
		util_validate()
		estimator=new("estimator", java=.jnew("com.portfolioeffect.quant.client.portfolio.Estimator",asset@java))
		return(estimator)
	}
	if(!exists('estimator')){
		stop('Could not create estimator object.')
	}
}

estimator_metric<-function(argList,estimator,...){
	util_validate(argList)
	data=list(...)
		resultTemp<-.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",toJSONpe(data))
		result<-getResult(resultTemp)
	return(result)
}


#price_ivRolling<-function(estimator,estimatorType=c("rv","tsrv","mrv","msrv","krv","jrmrv"),wLength=23400,...){
#	util_validate()
#	result<-switch(estimatorType[1],
#			rv=rvRolling(estimator,wLength=wLength),
#			tsrv=tsrvRolling(estimator,wLength=wLength,...),
#			mrv=mrvRolling(estimator,wLength=wLength),
#			msrv=msrvRolling(estimator,wLength=wLength,...),
#			krv=krvRolling(estimator,wLength=wLength,...),
#			jrmrv=jrmrvRolling(estimator,wLength=wLength))
#	return(result)
#}
price<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric='PRICE')
	return(result)
}

variance_rvRolling<-function(estimator,wLength=23400){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric='Rol_RV',wLength=as.integer(wLength))
	return(result)
}
	
variance_tsrvRolling<-function(estimator,K=2,wLength=23400){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="Rol_TSRV",numSubsamples = as.integer(K),wLength =as.integer(wLength))
	return(result)
}

variance_jrmrvRolling<-function(estimator,wLength=23400){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="Rol_JRMRV", wLength =as.integer(wLength))
	return(result)
}
	
variance_mrvRolling<-function(estimator,wLength=23400){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="Rol_MRV", wLength =as.integer(wLength))
	return(result)
}

variance_msrvRolling<-function(estimator,K=2,J=1,wLength=23400){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="Rol_MSRV",wLength = as.integer(wLength),  K = as.integer(K),  J = as.integer(J))
	return(result)
}

variance_krvRolling<-function(estimator,kernelName="ParzenKernel",bandwidth=1,wLength=23400){
	nameKernelJava<-NULL
	nameKernelJava<-switch(kernelName,
			BartlettKernel= paste("BARTLETT"),
			CubicKernel=paste("CUBIC"),	 
			EighthOrderKernel=paste("EIGHTH_ORDER"),	 
			EpanichnikovKernel=paste("EPANECHNIKOV"),	 
			FifthOrderKernel=paste("FIFTH_ORDER"),  
			ParzenKernel=paste("PARZEN"), 
			SecondOrderKernel=paste("SECOND_ORDER"), 
			SeventhOrderKernel=paste("SEVENTH_ORDER"), 
			SixthOrderKernel=paste("SIXTH_ORDER"), 
			TukeyHanningKernel=paste("TUKEY_HANNING"), 
			TukeyHanningModifiedKernel=paste("TUKEY_HANNING_MOD"),
			stop("Kernel not supported"))
	
#	if(bandwidth=="optimal"){
#		bandwidth<-computeKrvBandwidth(estimator,nameKernelJava)
#	}
	if(bandwidth<=0){
		stop("Incorrect value bandwidth")
	}
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="Rol_KRV",wLength = as.integer(wLength), kernelTypeName = nameKernelJava, bandwidth =as.integer(bandwidth))
	return(result)
}

#
# Noise variance
#

#price_nv<-function(estimator,estimatorType=c("rnv","acnv","urnv","uznv")){
#	util_validate()
#	result<-switch(estimatorType[1],
#			urnv=urnv(estimator),
#			acnv=acnv(estimator),
#			rnv=rnv(estimator),
#			uznv=uznv(estimator))
#
#	return(result)
#}

noise_urnv<-function(estimator)
{
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="URNV")
	return(result)
}

noise_acnv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="ACNV")
	return(result)
}

noise_rnv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="RNV")
	return(result)
}

noise_uznv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="UZNV")
	return(result)
}

#
# Integrated variance
#

#price_iv<-function(estimator,estimatorType=c("rv","tsrv","mrv","msrv","krv","jrmrv","uzrv"),...){
#	util_validate()
#	result<-switch(estimatorType[1],
#			rv=rv(estimator),
#			tsrv=tsrv(estimator,...),
#			mrv=mrv(estimator),
#			msrv=msrv(estimator),
#			krv=krv(estimator,...),
#			jrmrv=jrmrv(estimator),
#			uzrv=uzrv(estimator))
#
#	return(result)
#}

variance_rv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="RV")
	return(result)
}

variance_tsrv<-function(estimator,K=2){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="TSRV", numSubsamples = as.integer(K))
	return(result)
}

variance_mrv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="MRV")
	return(result)
}

variance_msrv<-function(estimator,K=2,J=1){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="MSRV",K = as.integer(K),J = as.integer(J))
	return(result)
}

variance_jrmrv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="JRMRV")
	return(result)
}

variance_krv<-function(estimator,kernelName="ParzenKernel",bandwidth=1){
	nameKernelJava<-NULL
	nameKernelJava<-switch(kernelName,
			BartlettKernel= paste("BARTLETT"),
			CubicKernel=paste("CUBIC"),	 
			EighthOrderKernel=paste("EIGHTH_ORDER"),	 
			EpanichnikovKernel=paste("EPANECHNIKOV"),	 
			FifthOrderKernel=paste("FIFTH_ORDER"),  
			ParzenKernel=paste("PARZEN"), 
			SecondOrderKernel=paste("SECOND_ORDER"), 
			SeventhOrderKernel=paste("SEVENTH_ORDER"), 
			SixthOrderKernel=paste("SIXTH_ORDER"), 
			TukeyHanningKernel=paste("TUKEY_HANNING"), 
			TukeyHanningModifiedKernel=paste("TUKEY_HANNING_MOD"),
			stop("Kernel not supported"))
	
	if(bandwidth<=0){
		stop("Incorrect value bandwidth")
	}
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="KRV",kernelTypeName = nameKernelJava,bandwidth = as.integer(bandwidth))
	return(result)
}

variance_uzrv<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="UZRV")
	return(result)
}

#
# Integrated quarticity
#

#price_iq<-function(estimator,estimatorType=c("rq","rqq","mrq","rtq","mtq")){
#	util_validate()
#	result<-switch(estimatorType[1],
#			rq=rq(estimator),
#			rqq=rqq(estimator),
#			mrq=mrq(estimator),
#			rtq=rtq(estimator),
#			mtq=mtq(estimator))
#
#	return(result)
#}

quarticity_rq<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="RQ")
	return(result)
}
quarticity_rqq<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="RQQ")
	return(result)
}
quarticity_mrq<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="MRQ")
	return(result)
}
quarticity_rtq<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="RTQ")
	return(result)
}
quarticity_mtq<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="MTQ")
	return(result)
}

#
# Noise-to-signal ratio
#

noise_nts<-function(estimator){
	result<-estimator_metric(argList=as.list(environment()),estimator=estimator,metric="NTS")
	return(result)}

estimator_availableSymbols<-function(estimator){
	result=.jcall(estimator@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getAllSymbolsList")
	result<-getResult(result)
	return(result)
}