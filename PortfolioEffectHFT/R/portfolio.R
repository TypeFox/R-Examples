setClass("portfolio",
		slots = c(java="jobjRef",optimization_info="ANY")
)

setMethod ("show" , "portfolio" ,
			function (object){
			util_validate()
			portfolio=portfolio_create(object)
			set<-portfolio_getSettings(portfolio)
			r<-NULL
			r<-rbind(r,c('Portfolio metrics mode',  	set$portfolioMetricsMode))
			r<-rbind(r,c('Window length ',				set$windowLength))
			r<-rbind(r,c('Time scale ',					set$timeScale))
			r<-rbind(r,c('Holding periods only ',		set$holdingPeriodsOnly))
			r<-rbind(r,c('Short sales mode',			set$shortSalesMode))
			r<-rbind(r,c('Price jumps model ',			set$jumpsModel))
			r<-rbind(r,c('Microstructure noise model ',	set$noiseModel))
			r<-rbind(r,c('Fractal price model ',	    set$fractalPriceModel))			
			r<-rbind(r,c('Portfolio factor model ',		set$factorModel))
			r<-rbind(r,c('Density model ',       		set$densityModel))
			r<-rbind(r,c('Drift term enabled ',        	set$driftTerm))
			r<-rbind(r,c('Results NA filter ',        	set$resultsNAFilter))
			r<-rbind(r,c('Results sampling interval ',	set$resultsSamplingInterval))
			r<-rbind(r,c('Input sampling interval ',	set$inputSamplingInterval))
			r<-rbind(r,c('Transaction cost per share ',	set$txnCostPerShare))
			r<-rbind(r,c('Transaction cost fixed ',     set$txnCostFixed))			
			
		colnames(r)<-array("",dim=2)
		rownames(r)<-array("",dim=NROW(r))
		cat (paste("PORTFOLIO SETTINGS",sep="    "))
			print(r, quote=FALSE)
			
			symbols<-portfolio_symbols(portfolio)
		
			.jcall(portfolio@java,returnSig="V", method="setParam","samplingInterval","last")
			
			portfolio_startBatch(portfolio)
			.jcall(portfolio@java,returnSig="V", method="createCallGroup",as.integer(7))
			.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",paste('{"metric":"PORTFOLIO_PROFIT"}'))
			.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",paste('{"metric":"PORTFOLIO_RETURN"}'))
			.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",paste('{"metric":"PORTFOLIO_VALUE"}'))
			
			portfolio_endBatch(portfolio)
			
			printMatrix1<-matrix(0,ncol=6,nrow=length(symbols))
#			j=1
#			for(symbol in symbols){
				printMatrix1[,1]<-position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_QUANTITY')[-1]
				printMatrix1[,2]<-round(position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_WEIGHT')[-1]*100,digits =3)
				printMatrix1[,3]<-round(position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_PROFIT')[-1],digits =3)
				printMatrix1[,4]<-round(position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_RETURN')[-1],digits =3)*100
				printMatrix1[,5]<-position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_VALUE')[-1]
				printMatrix1[,6]<-position_metric(argList=NULL,portfolio=portfolio,position='all',metric='POSITION_PRICE')[-1]
#				j<-j+1		
#			}
			dimnames(printMatrix1) = list( symbols, c("  Quantity","  Weight (in %)","  Profit","  Return (in %)","  Value","  Price"))
			symbols<-symbols[order(printMatrix1[,"  Weight (in %)"],decreasing=TRUE)]
			printMatrix1<-printMatrix1[symbols,]

			printMatrix2<-matrix(0,ncol=3,nrow=1)
			printMatrix2[1,1]<-round(portfolio_profit(portfolio)[2],digits =3)
			printMatrix2[1,2]<-round(portfolio_return(portfolio)[2],digits =3)*100
			printMatrix2[1,3]<-portfolio_value(portfolio)[2]
			dimnames(printMatrix2) = list( "Portfolio", c("  Profit","  Return (in %)","  Value"))
			cat (paste("\n\n","POSITION SUMMARY","\n",sep="    "))
			print(printMatrix1)
			cat (paste("\n\n","PORTFOLIO SUMMARY","\n",sep=""))
			print(printMatrix2)
})


util_summaryPlot<-function (x,y){	
util_summary(x)
}
#util_summaryPlotBW<-function (x,y){	
#	util_summary(x,bw=y)
#}
setMethod("plot" ,c(x="portfolio",y="missing"),util_summaryPlot)
#setMethod("plot" ,c(x="portfolio",y='logical'),util_summaryPlotBW)

portfolio_defaultSettings<-function(portfolio){
	portfolio_settings(portfolio,portfolioMetricsMode="portfolio",
				windowLength = "1d",
				holdingPeriodsOnly = FALSE,
				shortSalesMode = "lintner",
				jumpsModel = "moments",
				noiseModel = TRUE,
				fractalPriceModel=TRUE,
				factorModel = "sim",
				densityModel="GLD",
				driftTerm=FALSE,
				resultsNAFilter= TRUE,
				resultsSamplingInterval = "1s",
				inputSamplingInterval="1s",
				timeScale="1d",
				txnCostPerShare=0,
				txnCostFixed=0)
}
portfolio_settings<-function(
		portfolio,...
		){
	util_validate()
	isList=FALSE
	try({if(is.list(...)){isList=TRUE}},silent = TRUE)
	if(isList){change=(...)}else{change=list(...)}
	removeList<-NULL
	for(i in 1:length(change)){
		if(!(names(change[i]) %in% c( "portfolioMetricsMode","windowLength","holdingPeriodsOnly",     
								"shortSalesMode","jumpsModel","noiseModel", "fractalPriceModel",             
								"factorModel","resultsSamplingInterval","inputSamplingInterval",  
								"timeScale","driftTerm","txnCostPerShare",        
								"txnCostFixed","densityModel",'resultsNAFilter' ))){
				stopMessage('WRONG_SETTINGS_ARGUMENTS')
			}
		if(names(change[i]) %in% c( 'fractalPriceModel','holdingPeriodsOnly','noiseModel','resultsSamplingInterval','inputSamplingInterval','driftTerm','densityModel')){
		switch(names(change[i]),
				holdingPeriodsOnly = {
					change$isHoldingPeriodEnabled=if(change[[i]]){"true"}else{"false"}	
				},
				noiseModel={
					change$isNoiseModelEnabled=if(change[[i]]){"true"}else{"false"}
				},
				driftTerm = {
					change$isDriftEnabled=if(change[[i]]){"true"}else{"false"}
				},
				fractalPriceModel = {
					change$isFractalPriceModelEnabled=if(change[[i]]){"true"}else{"false"}
				},
				resultsSamplingInterval = {
					if(is.character(change[[i]])){
					change$samplingInterval=change[[i]]
				}else{
					change$samplingInterval=.jlong(change[[i]])
				}},
				inputSamplingInterval  = {
					change$priceSamplingInterval=change[[i]]
				},
				densityModel={
					change$densityApproxModel=change[[i]]
				}
				)
		removeList<-c(removeList,i)
}
}
change[removeList]<-NULL
for(i in 1:length(change)){
	if(names(change[i]) %in% 	c( 'resultsNAFilter')){
		.jcall(portfolio@java,returnSig="V", method="setNaNFiltered",as.logical(change[i]))
	}else{
	.jcall(portfolio@java,returnSig="V", method="setParam",names(change[i]),paste(change[[i]]))
}
}
}

portfolio_settingsRiskMetrics<-function(portfolio){
	portfolio_settings(portfolio,portfolioMetricsMode="price",
			windowLength = "1s",
			holdingPeriodsOnly = FALSE,
			shortSalesMode = "lintner",
			jumpsModel = "none",
			noiseModel = FALSE,
			fractalPriceModel=FALSE,
			factorModel = "sim",
			densityModel="NORMAL",
			driftTerm=FALSE,
			resultsSamplingInterval = "1d",
			inputSamplingInterval="1d",
			timeScale="1d",
			txnCostPerShare=0,
			txnCostFixed=0)
	.jcall(portfolio@java,returnSig="V", method="setParam","riskMethodology", "RiskMetrics")
}

portfolio_getSettings<-function(portfolio){

	temp<-list()
	temp$portfolioMetricsMode <- .jcall(portfolio@java,returnSig="S", method="getParam","portfolioMetricsMode")
	temp$windowLength<-.jcall(portfolio@java,returnSig="S", method="getParam","windowLength")
	temp$holdingPeriodsOnly<-as.logical(.jcall(portfolio@java,returnSig="S", method="getParam","isHoldingPeriodEnabled"))
	temp$shortSalesMode<-.jcall(portfolio@java,returnSig="S", method="getParam","shortSalesMode")
	temp$jumpsModel<-.jcall(portfolio@java,returnSig="S", method="getParam","jumpsModel")
	temp$noiseModel<-as.logical(.jcall(portfolio@java,returnSig="S", method="getParam","isNoiseModelEnabled"))
	temp$fractalPriceModel<-as.logical(.jcall(portfolio@java,returnSig="S", method="getParam","isNoiseModelEnabled"))
	temp$factorModel<-.jcall(portfolio@java,returnSig="S", method="getParam","factorModel")
	temp$resultsNAFilter<-.jcall(portfolio@java,returnSig="Z", method="isNaNFiltered")
	temp$resultsSamplingInterval<-.jcall(portfolio@java,returnSig="S", method="getParam","samplingInterval")
	temp$inputSamplingInterval<-.jcall(portfolio@java,returnSig="S", method="getParam","priceSamplingInterval")
	temp$timeScale<-.jcall(portfolio@java,returnSig="S", method="getParam","timeScale")
	temp$driftTerm<-as.logical(.jcall(portfolio@java,returnSig="S", method="getParam","isDriftEnabled"))
	temp$txnCostPerShare<-.jcall(portfolio@java,returnSig="S", method="getParam","txnCostPerShare")
	temp$txnCostFixed<-.jcall(portfolio@java,returnSig="S", method="getParam","txnCostFixed")	
	temp$densityModel<-.jcall(portfolio@java,returnSig="S", method="getParam","densityApproxModel")
	temp
}

portfolio_create<-function(index=NULL,fromTime=NULL,toTime=NULL,priceDataIx=NULL){
	
	if(is.null(index) & is.null(fromTime)& is.null(toTime)& is.null(priceDataIx)){
		stop('No arguments provided, please check required arguments list.')
	}
	if(((class(index)=="matrix")&(is.null(priceDataIx)))||((class(priceDataIx)=="matrix")&(is.null(index)))){
		if((class(index)=="matrix")&(is.null(priceDataIx))){
			priceDataIx=index
		}
		util_validate()
		clientConnection=getOption('clientConnection')
		portfolio=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",clientConnection),optimization_info=NULL)
		result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addIndex", as.double(priceDataIx[,2]),.jlong(priceDataIx[,1]))
		util_checkErrors(result)
		portfolio_defaultSettings(portfolio)
		return(portfolio)
	}
	if(((class(fromTime)=="character")&(class(toTime)=="character")&(is.null(priceDataIx)))||((class(fromTime)=="character")&(class(index)=="character")&(is.null(toTime))&(is.null(priceDataIx)))){
		if(is.null(index)){
			index="SPY"
		}
		if(((class(fromTime)=="character")&(class(index)=="character")&(is.null(toTime))&(is.null(priceDataIx)))){
			toTime=fromTime
			fromTime=index
			index="SPY"
		}
		util_validate()
		clientConnection=getOption('clientConnection')
		portfolio=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",clientConnection),optimization_info=NULL)
		result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setFromTime", fromTime)
		util_checkErrors(result)
		result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setToTime", toTime)
		util_checkErrors(result)
		result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addIndex", index)
		util_checkErrors(result)
		portfolio_defaultSettings(portfolio)
		return(portfolio)
	}
	if((class(index)=="portfolio")&(is.null(fromTime))&(is.null(toTime))&(is.null(priceDataIx))){
		util_validate()
		portfolio=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",index@java),optimization_info=NULL)
		return(portfolio)
	}
	if(!exists('portfolio')){
		stop('Could not create portfolio object.')
	}
}

portfolio_addPosition<-function(portfolio,symbol,quantity,time,priceData){
}

setMethod("portfolio_addPosition" ,c(portfolio="portfolio",symbol="character",quantity="ANY",time="missing",priceData="missing"),function(
				portfolio,symbol,quantity){
			util_validate()
			result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addPosition",symbol, as.integer(quantity))
			util_checkErrors(result)
		})

setMethod("portfolio_addPosition" ,c(portfolio="portfolio",symbol="character",quantity="ANY",time="missing",priceData="matrix"),function(
				portfolio,symbol,quantity,priceData){
			util_validate()
			result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addPosition",symbol, as.double(priceData[,2]), as.integer(quantity),.jlong(priceData[,1]))
			util_checkErrors(result)
		})

setMethod("portfolio_addPosition" ,c(portfolio="portfolio",symbol="character",quantity="ANY",time="ANY",priceData="missing"),function(
				portfolio,symbol,quantity,time){
			util_validate()
			if(!is.numeric(time)){
				time<-util_dateToPOSIXTime(time)
			}
			result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addPosition",symbol, as.integer(quantity),.jlong(time))
			util_checkErrors(result)
		})

setMethod("portfolio_addPosition" ,c(portfolio="portfolio",symbol="character",quantity="ANY",time="ANY",priceData="matrix"),function(
				portfolio,symbol,quantity,time,priceData){
			util_validate()
	if(!is.numeric(time)){
		time<-util_dateToPOSIXTime(time)
	}
			result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="addPosition",symbol, as.double(priceData[,2]),.jlong(priceData[,1]), as.integer(quantity),.jlong(time))
			util_checkErrors(result)
		})


portfolio_removePosition<-function(portfolio,symbol){
	util_validate(as.list(environment()))
	if(!is(portfolio,'portfolio')){
		stopMessage("NOT_PORTFOLIO_CLASS")
	}
	if(!is(symbol,"character")){
		stopMessage("SYMBOL_NOT_CHARACTER_CLASS")
	}
	.jcall(portfolio@java,returnSig="V", method="removePositionQuantity",symbol)
	.jcall(portfolio@java,returnSig="V", method="removePositionPrice",symbol)
}

position_setQuantity<-function(portfolio,symbol,quantity){
	util_validate()
	result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setPositionQuantity",symbol, as.integer(quantity))  
	util_checkErrors(result)
}

position_quantity<-function(portfolio,symbol){
	util_validate(as.list(environment()))
	result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getPositionQuantity",symbol)  
	util_checkErrors(result)
	result<-getResultValuesDoubleArrayWithTime(result)
	return(result)}

portfolio_symbols<-function(portfolio){
	util_validate(as.list(environment()))
	result<-.jcall(portfolio@java,returnSig="[S", method="getSymbols")  
	return(result)}


toJSONpe<-function(x){
	names=names(x)
	result='{'
	t=0
	for(name in names){
		t=t+1
		if(t<=1){
			result=paste(result,'"',name,'"',':','"',x[[name]],'"',sep="")
		}else{
			result=paste(result,', "',name,'"',':','"',x[[name]],'"',sep="")
		}
	}
	result=paste(result,'}',sep="")
	return(result)
}

position_metric<-function(argList,portfolio,...){
	util_validate(argList)
	data=list(...)
	result=NULL
	if((data$position=='all')&&(!is.null(data$position))){
		set<-portfolio_getSettings(portfolio)
		if(set$resultsSamplingInterval=="last"){
			data[["position"]]<-NULL
	data$value=data$metric
	data$metric="POSITION_MATRIX"
			resultTemp<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",toJSONpe(data))
			result<-getResult(resultTemp)
	}else{
		symbols=portfolio_symbols(portfolio)
		portfolio_startBatch(portfolio)
		for(symbol in symbols){
			data$position=symbol
			.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",toJSONpe(data))
		}
		portfolio_endBatch(portfolio)
		for(symbol in symbols){
			data$position=symbol
			resultTemp<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",toJSONpe(data))
		util_checkErrors(resultTemp)
		result<-cbind(result,getResultValuesDoubleArray(resultTemp))
	}
	result<-cbind(getTimeMilliSec(resultTemp),result)
	colnames(result)<-c("Time",symbols)
}
	}else{
		resultTemp<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getMetric",toJSONpe(data))
		result<-getResult(resultTemp)

	}
	return(result)
}

position_maxDrawdown<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_MAX_DRAWDOWN')
	return(result)}

position_calmarRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_CALMAR_RATIO')
	return(result)}

position_profit<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_PROFIT')
	return(result)}

position_price<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_PRICE')
	return(result)}

position_weight<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_WEIGHT')
	return(result)}

position_value<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_VALUE')
	return(result)}

position_beta<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_BETA')
	return(result)}

position_alpha<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_ALPHA')
	return(result)}

position_return<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_RETURN')
	return(result)}

position_expectedReturn<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_EXPECTED_RETURN')
	return(result)
}

position_returnAutocovariance<-function(portfolio,symbol,lag=10){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_RETURN_AUTOCOVARIANCE',lag=lag)
	return(result)
}

position_variance<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_VARIANCE')
	return(result)
}

position_CVaR<-function(portfolio,symbol,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_CVAR',confidenceInterval=as.double(confidenceInterval))
	return(result)}

position_VaR<-function(portfolio,symbol,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_VAR',confidenceInterval=as.double(confidenceInterval))
	return(result)}

position_modifiedSharpeRatio<-function(portfolio,symbol,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_SHARPE_RATIO_MOD',confidenceInterval=as.double(confidenceInterval))
	return(result)}

position_starrRatio<-function(portfolio,symbol,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_STARR_RATIO',confidenceInterval=as.double(confidenceInterval))
	return(result)}

position_sharpeRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_SHARPE_RATIO')
	return(result)}

position_treynorRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_TREYNOR_RATIO')
	return(result)}

position_skewness<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_SKEWNESS')
	return(result)}

position_kurtosis<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_KURTOSIS')
	return(result)}

position_informationRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_INFORMATION_RATIO')
	return(result)
}

position_jensensAlpha<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_ALPHA_JENSEN')
	return(result)
}

position_covariance<-function(portfolio,symbol1,symbol2){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='POSITION_COVARIANCE',positionA=symbol1,positionB=symbol2)
	return(result)
}

position_correlation<-function(portfolio,symbol1,symbol2){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='POSITION_CORRELATION',positionA=symbol1,positionB=symbol2)
	return(result)
}

position_omegaRatio<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_OMEGA_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}


position_rachevRatio<-function(portfolio,symbol,confidenceIntervalA=0.95,confidenceIntervalB=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_RACHEV_RATIO',confidenceIntervalAlpha=as.double(confidenceIntervalA),confidenceIntervalBeta=as.double(confidenceIntervalB))
	return(result)}

position_gainVariance<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_GAIN_VARIANCE')
	return(result)}

position_lossVariance<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_LOSS_VARIANCE')
	return(result)}

position_downsideVariance<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_DOWNSIDE_VARIANCE',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_upsideVariance<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_UPSIDE_VARIANCE',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_expectedDownsideReturn<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_EXPECTED_DOWNSIDE_THRESHOLD_RETURN',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_expectedUpsideReturn<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_EXPECTED_UPSIDE_THRESHOLD_RETURN',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_sortinoRatio<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_SORTINO_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_upsideDownsideVarianceRatio<-function(portfolio,symbol,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_UPSIDE_DOWNSIDE_VARIANCE_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}

position_gainLossVarianceRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_GAIN_LOSS_VARIANCE_RATIO')
	return(result)}

position_downCaptureRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_DOWN_CAPTURE_RATIO')
	return(result)}

position_upCaptureRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_UP_CAPTURE_RATIO')
	return(result)}

position_downNumberRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_DOWN_NUMBER_RATIO')
	return(result)}

position_upNumberRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_UP_NUMBER_RATIO')
	return(result)}

position_downPercentageRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_DOWN_PERCENTAGE_RATIO')
	return(result)}

position_hurstExponent<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_HURST_EXPONENT')
	return(result)}

position_fractalDimension<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_FRACTAL_DIMENSION')
	return(result)}

position_moment<-function(portfolio,symbol,order){
	util_validate()
	if(order=="all"){
		order=c(3,4)
	}
	totalResult<-NULL
	for(i in order){
		result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric=paste('POSITION_MOMENT',i,sep=""))
		if(is.null(totalResult)){
			totalResult<-result
		}else{
			totalResult=cbind(totalResult,result[,2])
		}
	}
	return(totalResult)}


position_cumulant<-function(portfolio,symbol,order){
	util_validate()
	if(order=="all"){
		order=c(3,4)
	}
	totalResult<-NULL
	for(i in order){
		result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric=paste('POSITION_CUMULANT',i,sep=""))
		if(is.null(totalResult)){
			totalResult<-result
		}else{
			totalResult=cbind(totalResult,result[,2])
		}
	}
	return(totalResult)}

position_upPercentageRatio<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_UP_PERCENTAGE_RATIO')
	return(result)}

position_covarianceMatrix<-function(portfolio){
	util_validate()
	portfolioTemp=portfolio_create(portfolio)
	symbols<-portfolio_symbols(portfolioTemp)
	
    .jcall(portfolioTemp@java,returnSig="V", method="setSamplingInterval","last")
		
	n<-length(symbols)
	resultMatrix<-matrix(1,nrow=n,ncol=n)
	.jcall(portfolioTemp@java,returnSig="V", method="createCallGroup",as.integer(n+((n*n-n)/2)))
	for(i in 2:n){
		for(j in 1:(i-1)){
			resultMatrix[i,j]<-position_metric(argList=as.list(environment()),portfolio=portfolioTemp,metric='POSITION_COVARIANCE',positionA=symbols[i],positionB=symbols[j])[2]
			resultMatrix[j,i]<-resultMatrix[i,j]
		}
	}
	for(i in 1:n){

		resultMatrix[i,i]<-position_metric(argList=as.list(environment()),portfolio=portfolioTemp,position=symbols[i],metric='POSITION_VARIANCE')[2]
	}
	dimnames(resultMatrix) = list(symbols,symbols)
	return(resultMatrix)
}

position_correlationMatrix<-function(portfolio){
	util_validate()
	symbols<-portfolio_symbols(portfolio)
	n<-length(symbols)
	
	portfolioTemp=portfolio_create(portfolio)
    .jcall(portfolioTemp@java,returnSig="V", method="setSamplingInterval","last")
		
	resultMatrix<-matrix(1,nrow=n,ncol=n)
	.jcall(portfolioTemp@java,returnSig="V", method="createCallGroup",as.integer(((n*n-n)/2)))
	for(i in 2:n){
		for(j in 1:(i-1)){
			resultMatrix[i,j]<-position_metric(argList=as.list(environment()),portfolio=portfolioTemp,metric='POSITION_CORRELATION',positionA=symbols[i],positionB=symbols[j])[2]
			resultMatrix[j,i]<-resultMatrix[i,j]
		}
	}
	dimnames(resultMatrix) = list(symbols,symbols)
	return(resultMatrix)
}

position_txnCosts<-function(portfolio,symbol){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,position=symbol,metric='POSITION_TRANSACTION_COSTS_SIZE')
	return(result)}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~ Portfolio Methods ~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

portfolio_value<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_VALUE')
	return(result)}

portfolio_return<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_RETURN')
	return(result)}

portfolio_expectedReturn<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_EXPECTED_RETURN')
	return(result)}

portfolio_profit<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_PROFIT')
	return(result)}

portfolio_beta<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_BETA')
	return(result)}

portfolio_alpha<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_ALPHA')
	return(result)}

portfolio_variance<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_VARIANCE')
	return(result)}

portfolio_maxDrawdown<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_MAX_DRAWDOWN')
	return(result)}

portfolio_calmarRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_CALMAR_RATIO') 
	return(result)}

portfolio_VaR<-function(portfolio,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_VAR',confidenceInterval=as.double(confidenceInterval))
	return(result)}

portfolio_CVaR<-function(portfolio,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_CVAR',confidenceInterval=as.double(confidenceInterval))
	return(result)}

portfolio_modifiedSharpeRatio<-function(portfolio,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_SHARPE_RATIO_MOD',confidenceInterval=as.double(confidenceInterval))
	return(result)}

portfolio_starrRatio<-function(portfolio,confidenceInterval=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_STARR_RATIO',confidenceInterval=as.double(confidenceInterval))
	return(result)}

portfolio_sharpeRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_SHARPE_RATIO')
	return(result)}

portfolio_treynorRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_TREYNOR_RATIO')
	return(result)}

portfolio_skewness<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_SKEWNESS')
	return(result)}

portfolio_kurtosis<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_KURTOSIS')
	return(result)}

portfolio_informationRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_INFORMATION_RATIO')
	return(result)}

portfolio_jensensAlpha<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_ALPHA_JENSEN')
	return(result)}

portfolio_omegaRatio<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_OMEGA_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_rachevRatio<-function(portfolio,confidenceIntervalA=0.95,confidenceIntervalB=0.95){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_RACHEV_RATIO',confidenceIntervalAlpha=as.double(confidenceIntervalA),confidenceIntervalBeta=as.double(confidenceIntervalB))
	return(result)}

portfolio_gainVariance<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_GAIN_VARIANCE')
	return(result)}

portfolio_lossVariance<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_LOSS_VARIANCE')
	return(result)}

portfolio_downsideVariance<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_DOWNSIDE_VARIANCE',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_upsideVariance<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_UPSIDE_VARIANCE',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_expectedDownsideReturn<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_EXPECTED_DOWNSIDE_THRESHOLD_RETURN',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_expectedUpsideReturn<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_EXPECTED_UPSIDE_THRESHOLD_RETURN',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_hurstExponent<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_HURST_EXPONENT')
	return(result)}

portfolio_fractalDimension<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_FRACTAL_DIMENSION')
	return(result)}

portfolio_txnCosts<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_TRANSACTION_COSTS_SIZE')
	return(result)}

portfolio_moment<-function(portfolio,order){
	if(order=="all"){
		order=c(3,4)
	}
	totalResult<-NULL
	for(i in order){
		result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric=paste("PORTFOLIO_MOMENT",i,sep=""))
		if(is.null(totalResult)){
			totalResult<-result
		}else{
			totalResult=cbind(totalResult,result[,2])
		}
	}
	return(totalResult)}

portfolio_cumulant<-function(portfolio,order){
	if(order=="all"){
		order=c(3,4)
	}
	totalResult<-NULL
	for(i in order){
		result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric=paste("PORTFOLIO_CUMULANT",i,sep=""))
		if(is.null(totalResult)){
			totalResult<-result
		}else{
			totalResult=cbind(totalResult,result[,2])
		}
	}
	return(totalResult)}

portfolio_sortinoRatio<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_SORTINO_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_upsideDownsideVarianceRatio<-function(portfolio,thresholdReturn){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_UPSIDE_DOWNSIDE_VARIANCE_RATIO',thresholdReturn=as.double(thresholdReturn))
	return(result)}

portfolio_gainLossVarianceRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_GAIN_LOSS_VARIANCE_RATIO')
	return(result)}

portfolio_downCaptureRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_DOWN_CAPTURE_RATIO')
	return(result)}

portfolio_upCaptureRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_UP_CAPTURE_RATIO')
	return(result)}

portfolio_downNumberRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_DOWN_NUMBER_RATIO')
	return(result)}

portfolio_upNumberRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_UP_NUMBER_RATIO')
	return(result)}

portfolio_downPercentageRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_DOWN_PERCENTAGE_RATIO')
	return(result)}

portfolio_upPercentageRatio<-function(portfolio){
	result<-position_metric(argList=as.list(environment()),portfolio=portfolio,metric='PORTFOLIO_UP_PERCENTAGE_RATIO')
	return(result)}

portfolio_pdf<-function(portfolio,pValueLeft,pValueRight,nPoints,addNormalDensity=FALSE){
	
	portfolioTemp=portfolio_create(portfolio)
	set<-portfolio_getSettings(portfolioTemp)
    .jcall(portfolioTemp@java,returnSig="V", method="setSamplingInterval","last")
		
	z<-.jcall(portfolioTemp@java, returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getPDF",
			as.double(pValueLeft),as.double(pValueRight),as.integer(nPoints))
	
	result<-list(pdf=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "pdf" ,simplify=TRUE),
			value=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "x", simplify=TRUE),
			time=.jcall(z,returnSig="[J",method="getLongArray", "time"))
	
	if(addNormalDensity){
		GaussianMoments=FALSE
		if(set$densityModel!="NORMAL"){
			GaussianMoments<-TRUE
			portfolio_settings(portfolioTemp,densityModel="NORMAL")
		}
		z<-.jcall(portfolioTemp@java,
				returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getPDF",
				as.double(pValueLeft),as.double(pValueRight),as.integer(nPoints))
		
		result$pdfNormal=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "pdf", simplify=TRUE)
		result$valueNormal=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "x", simplify=TRUE)
	}
	return(result)
}

position_pdf<-function(portfolio,symbol,pValueLeft,pValueRight,nPoints,addNormalDensity=FALSE){
	
	portfolioTemp=portfolio_create(portfolio)
	set<-portfolio_getSettings(portfolioTemp)
    .jcall(portfolioTemp@java,returnSig="V", method="setSamplingInterval","last")
    
	z<-.jcall(portfolioTemp@java,
			returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getPDF",
			as.double(pValueLeft),as.double(pValueRight),as.integer(nPoints), symbol)
	
	result<-list(pdf=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "pdf", simplify=TRUE),
			value=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "x", simplify=TRUE),
			time=.jcall(z,returnSig="[J",method="getLongArray", "time"))
	
	if(addNormalDensity){
		GaussianMoments=FALSE
		if(set$densityModel!="NORMAL"){
			GaussianMoments<-TRUE
			portfolio_settings(portfolioTemp,densityModel="NORMAL")
		}
		z<-.jcall(portfolioTemp@java,
				returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getPDF",
				as.double(pValueLeft),as.double(pValueRight),as.integer(nPoints), symbol)
		
		result$pdfNormal=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "pdf", simplify=TRUE)
		result$valueNormal=.jcall(z,returnSig="[[D",method="getDoubleMatrix", "x", simplify=TRUE)
	}
	
	return(result)
}

position_returnJumpSize<-function(portfolio,symbol){
	portfolioTemp=portfolio_create(portfolio)
	portfolio_settings(portfolioTemp,jumpsModel='none')
    time=position_price(portfolioTemp,symbol)[,1]
	priceNoJumpsFilter=position_price(portfolioTemp,symbol)[,2]
	portfolio_settings(portfolioTemp,jumpsModel='all')
	priceJumpsFilter=position_price(portfolioTemp,symbol)[,2]
	jumps=log(priceNoJumpsFilter)-log(priceJumpsFilter)
	cbind(time,jumps)
}

portfolio_startBatch<-function(portfolio){
	.jcall(portfolio@java,returnSig="V", method="startBatch")
}

portfolio_endBatch<-function(portfolio){
	result=.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="finishBatch")
	util_checkErrors(result)
}


portfolio_availableSymbols<-function(portfolio){
	result=.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="getAllSymbolsList")
	result<-getResult(result)
	return(result)
}