setClass("optimizer",
		slots = c(java="jobjRef",
		portfolio="jobjRef",
		errorInDecimalPoints='numeric',
		globalOptimumProbability='numeric',
		constraintMerticFunctions='ANY',
		constraintTypeFunctions='ANY',
		functions='ANY',
		constraintConfidenceIntervalFunctions='ANY',
		constraintSymbolsFunctions='ANY',
		constraintMerticSimple='ANY',
		constraintTypeSimple='ANY',
		constraintValueSimple='ANY',
		constraintConfidenceInterval='ANY',
		constraintSymbols='ANY',
		portfolioValue='ANY',
		goal="character",
		direction="character",
		confidenceInterval='numeric',
		forecastedValueLists='ANY',
		forecastTimeStep='ANY',
		forecastType="character",
		forecastExponentialWindow="character",
		forecastPortfolioWindow="character"
))



optimization_run<-function(optimizer){
	util_validate(as.list(environment()))
	
	## --- BEGIN forecasting code -- ##
#	forecastingPortfolio=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",optimizer@portfolio),optimization_info=NULL)
	settings = portfolio_getSettings(new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",optimizer@portfolio),optimization_info=NULL))
#	.jcall(forecastingPortfolio@java,returnSig="V", method="setParam",'isRebalancingHistoryEnabled',"false")
#	.jcall(forecastingPortfolio@java,returnSig="V", method="setParam",'windowLength',optimizer@windowLength)
	

	if(settings$portfolioMetricsMode == 'portfolio') 
	{
		forecastedValues=.jnew("com.portfolioeffect.quant.client.portfolio.optimizer.ForecastedValues",optimizer@portfolio)
		
		result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setForecastTimeStep",optimizer@forecastTimeStep)
		util_checkErrors(result)
#		result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="makeSimpleCumulantsForecast",forecastingPortfolio@java,optimizer@forecastType)
#		util_checkErrors(result)
		.jcall(optimizer@java,returnSig="V",method="setForecasterType",optimizer@forecastType)
		.jcall(optimizer@java,returnSig="V",method="setForecastExpWindow",optimizer@forecastExponentialWindow)
		.jcall(optimizer@java,returnSig="V",method="setForecastPortfolioWindow",optimizer@forecastPortfolioWindow)
		if(!is.null(optimizer@forecastedValueLists))
		{
			for(forecastedValueList in optimizer@forecastedValueLists){	
				switch(forecastedValueList$metricType,
			ExpReturn={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedExpReturn",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Beta={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedBeta",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Variance={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedVariance",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Skewness={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedSkewness",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Kurtosis={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedKurtosis",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Cumulant1={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedCumulant1",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Cumulant2={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedCumulant2",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Cumulant3={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedCumulant3",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		Cumulant4={result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="setSymbolForecastedCumulant4",forecastedValueList$symbol,as.double(forecastedValueList$value),.jlong(util_dateToPOSIXTime(forecastedValueList$time)))
			util_checkErrors(result)},
		stop("Incorrect forecast metric type"))
				# TODO finish with the list
				# result<-.jcall(forecastedValues,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="makeSimpleCumulantsForecast")
				# util_checkErrors(result)	
			}
		}
		.jcall(optimizer@java,returnSig="V",method="setForecastedValue",forecastedValues)
	}
	## --- END forecasting code -- ##
	
	if(is.null(optimizer@functions)){
		result<-.jcall(optimizer@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getOptimizedPortfolio")
		portfolio<-new("portfolio", java=getResult(result),optimization_info=NULL)
		portfolio@optimization_info=c(.jcall(result,returnSig="S", method="getInfoParam","nFunctionExecute"),.jcall(result,returnSig="S", method="getInfoParam","nGlobalStart"),.jcall(result,returnSig="S", method="getInfoParam","nLocalSolution"),.jcall(result,returnSig="S", method="getInfoParam","nOptimizations"),.jcall(result,returnSig="S", method="getInfoParam","nConstraintSatisfied"))
	}else{
		portfolio<-util_optimizationFunction(optimizer)
	}
	return(portfolio)
}

optimization_forecast<-function(optimizer,metricType,symbol,value,time){
	util_validate(as.list(environment()))
		optimizer@forecastedValueLists=c(optimizer@forecastedValueLists,list(list(metricType=metricType,symbol=symbol,value=value,time=time)))
		return(optimizer)
}

optimization_constraint<-function(optimizer,
		constraintType,
		constraintMertic,
		constraintValue,confidenceInterval=NULL,symbols=NULL){
	constraintTypeFinal = switch(constraintType[1],
			"=" = "equals",
			">=" = "greaterOrEquals",
			"<=" = "lessOrEquals",
			"equals" = "equals",
			"greaterOrEquals" = "greaterOrEquals",
			"lessOrEquals" = "lessOrEquals",
			stop("Incorrect position constraint type")
	)
	if(class(constraintValue)=="function"){
		function_names=names(formals(constraintValue));
#			if('portfolio' %in% function_names){
		optimizer@constraintMerticFunctions=c(optimizer@constraintMerticFunctions,list(constraintMertic))
		optimizer@constraintTypeFunctions=c(optimizer@constraintTypeFunctions,list(constraintTypeFinal))
		optimizer@functions=c(optimizer@functions,list(constraintValue))
		optimizer@constraintConfidenceIntervalFunctions=c(optimizer@constraintConfidenceIntervalFunctions,list(confidenceInterval))
		optimizer@constraintSymbolsFunctions=c(optimizer@constraintSymbolsFunctions,list(symbols))
#		}else{if('time' %in% function_names){
#				CV=NULL
#				for(time in optimizer@rebalancingTimes){
#					CV=c(CV,constraintValue(time))
#				}
#					.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(CV),.jlong(util_dateToPOSIXTime(optimizer@rebalancingTimes)))
#					optimizer@constraintValueSimple=c(optimizer@constraintValueSimple,list(cbind(as.double(CV),util_dateToPOSIXTime(optimizer@rebalancingTimes))))
#				}else{
#					CV=NULL
#					for(time in optimizer@rebalancingTimes){
#						CV=c(CV,constraintValue())
#					}
#				.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(CV),.jlong(util_dateToPOSIXTime(optimizer@rebalancingTimes)))	
#				optimizer@constraintValueSimple=c(optimizer@constraintValueSimple,list(cbind(as.double(CV),util_dateToPOSIXTime(optimizer@rebalancingTimes))))
#			}
#			optimizer@constraintMerticSimple=c(optimizer@constraintMerticSimple,list(constraintMertic))
#			optimizer@constraintTypeSimple=c(optimizer@constraintTypeSimple,list(constraintTypeFinal))
#		}
	}else{
		if(is.null(ncol(constraintValue))){
			switch(constraintMertic,
					BETA=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue)),
					POSITION_WEIGHT=
							.jcall(optimizer@java,returnSig="V",method="addPositionConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),if(!is.null(symbols)){symbols}),
					RETURN=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue)),
					EXPECTED_RETURN=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue)),
					VARIANCE=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue)),
					SHARPE_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue)),
					MODIFIED_SHARPE_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),as.double(confidenceInterval)),
					STARR_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),as.double(confidenceInterval)),
					VAR=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),as.double(confidenceInterval)),
					CVAR=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),as.double(confidenceInterval)),
					POSITIONS_SUM_ABS_WEIGHT=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue),symbols))
			
		}else{
			switch(constraintMertic,
					BETA=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1]))),
					POSITION_WEIGHT=
							.jcall(optimizer@java,returnSig="V",method="addPositionConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),if(!is.null(symbols)){symbols}),
					RETURN=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1]))),
					EXPECTED_RETURN=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1]))),
					VARIANCE=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1]))),
					SHARPE_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1]))),
					MODIFIED_SHARPE_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),as.double(confidenceInterval)),
					STARR_RATIO=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),as.double(confidenceInterval)),
					VAR=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),as.double(confidenceInterval)),
					CVAR=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),as.double(confidenceInterval)),
					POSITIONS_SUM_ABS_WEIGHT=
							.jcall(optimizer@java,returnSig="V",method="addPortfolioConstraint",constraintMertic,constraintTypeFinal,as.double(constraintValue[,2]),.jlong(util_dateToPOSIXTime(constraintValue[,1])),symbols))
		}
		optimizer@constraintValueSimple=c(optimizer@constraintValueSimple,list(constraintValue))
		optimizer@constraintMerticSimple=c(optimizer@constraintMerticSimple,list(constraintMertic))
		optimizer@constraintTypeSimple=c(optimizer@constraintTypeSimple,list(constraintTypeFinal))
		optimizer@constraintConfidenceInterval=c(optimizer@constraintConfidenceInterval,list(confidenceInterval))
		optimizer@constraintSymbols=c(optimizer@constraintSymbols,list(symbols))
	}
	return(optimizer)
}


optimization_goal<-function(portfolio,
		goal=c("EquiWeight",
				"ContraintsOnly",
				"Variance",
				"VaR",
				"CVaR",
				"ExpectedReturn",
				"Return",
				"SharpeRatio",
				"ModifiedSharpeRatio",
				"StarrRatio"
		),
		direction=c("minimize", 
				"maximize"),
		confidenceInterval=0.95,
		forecastPortfolioWindow='1m',
		forecastTimeStep='1m',
		forecastType=c("exp_smoothing", "simple"),
		forecastExponentialWindow='5m',
		errorInDecimalPoints=1e-12,
		globalOptimumProbability=0.99)
{
	util_validate(as.list(environment()))
if((globalOptimumProbability>=1)|(globalOptimumProbability<=0)){
	stop("globalOptimizationProbability is bounded by (0, 1) excluding interval endpoints")
}
	if(!(direction[1] %in% c("minimize", "maximize"))) {
		stop("Direction not specified")
	}
	if(.jcall(portfolio@java,returnSig="S", method="getParam","portfolioMetricsMode")=='portfolio'){
		path<-"com.portfolioeffect.quant.client.portfolio.optimizer.StrategyOptimizer"
	}else{
		path<-"com.portfolioeffect.quant.client.portfolio.optimizer.PortfolioOptimizer"
	}
	
		optimizer<-new("optimizer", java=.jnew(path,portfolio@java,as.double(errorInDecimalPoints),as.double(globalOptimumProbability)),
				portfolio=portfolio@java,errorInDecimalPoints=as.double(errorInDecimalPoints),globalOptimumProbability=as.double(globalOptimumProbability),constraintMerticFunctions=NULL,constraintTypeFunctions=NULL,functions=NULL,
				portfolioValue=NULL, goal=goal, direction=direction, confidenceInterval=confidenceInterval,
				constraintMerticSimple=NULL,	constraintTypeSimple=NULL,	constraintValueSimple=NULL,	constraintConfidenceInterval=NULL,
				constraintSymbols=NULL,	constraintConfidenceIntervalFunctions=NULL,constraintSymbolsFunctions=NULL,		forecastedValueLists=NULL,
				forecastTimeStep='1m',forecastType=forecastType[1],forecastExponentialWindow=forecastExponentialWindow,forecastPortfolioWindow=forecastPortfolioWindow)

	
	switch(goal[1], 
			Variance=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","VARIANCE",direction[1]),
			VaR=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","VAR", direction[1], as.double(confidenceInterval)),
			CVaR=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","CVAR", direction[1], as.double(confidenceInterval)),
			ExpectedReturn=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","EXPECTED_RETURN", direction[1]),
			Return=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","RETURN", direction[1]),
			SharpeRatio=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","SHARPE_RATIO", direction[1]),
			ModifiedSharpeRatio=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","SHARPE_RATIO", direction[1], as.double(confidenceInterval)),
			StarrRatio=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","STARR_RATIO", direction[1], as.double(confidenceInterval)),
			EquiWeight=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","NONE", direction[1]),
			ContraintsOnly=.jcall(optimizer@java,returnSig="V",method="setOptimizationGoal","ZERO", direction[1]),
			
			stop("Optimization goal not supported"))
	
	if(is.numeric(forecastTimeStep)){
		forecastTimeStep=paste(forecastTimeStep,'s',sep="")
	}
	optimizer@forecastTimeStep=forecastTimeStep
	return(optimizer)
}

optimization_forecasted_values<-function(optimizer,symbol,forecastedMetric,forecastedTime,forecastedValue){
	
}

optimization_constraint_portfolioValue<-function(
		optimizer,
		constraintValue){
	util_validate(as.list(environment()))
	.jcall(optimizer@java,returnSig="V",method="setPortfolioValue",as.integer(constraintValue));
	optimizer@portfolioValue=as.integer(constraintValue);
	return(optimizer)
}

optimization_constraint_beta<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	constraintMertic="BETA"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
return(optimizer)
}

optimization_constraint_allWeights<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	
	constraintMertic="POSITION_WEIGHT"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
	return(optimizer)
}

optimization_constraint_return<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	
	constraintMertic="RETURN"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
	return(optimizer)
}
optimization_constraint_expectedReturn<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	
	constraintMertic="EXPECTED_RETURN"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
	return(optimizer)
}
optimization_constraint_variance<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	
	constraintMertic="VARIANCE"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
	return(optimizer)
}

optimization_constraint_sharpeRatio<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue){
	util_validate(as.list(environment()))
	
	constraintMertic="SHARPE_RATIO"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue)
	return(optimizer)
}
optimization_constraint_modifiedSharpeRatio<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		confidenceInterval=0.95){
	util_validate(as.list(environment()))
	
	constraintMertic="MODIFIED_SHARPE_RATIO"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,confidenceInterval=confidenceInterval)
	return(optimizer)
}
optimization_constraint_starrRatio<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		confidenceInterval=0.95){
	util_validate(as.list(environment()))
	
	constraintMertic="STARR_RATIO"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,confidenceInterval=confidenceInterval)
	return(optimizer)
}
optimization_constraint_VaR<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		confidenceInterval=0.95){
	util_validate(as.list(environment()))
	
	constraintMertic="VAR"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,confidenceInterval=confidenceInterval)
	return(optimizer)
}
optimization_constraint_CVaR<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		confidenceInterval=0.95){
	util_validate(as.list(environment()))
	
	constraintMertic="CVAR"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,confidenceInterval=confidenceInterval)
	return(optimizer)
}


optimization_constraint_sumOfAbsWeights<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		symbols){
	util_validate(as.list(environment()))
	
	constraintMertic="POSITIONS_SUM_ABS_WEIGHT"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,symbols=symbols)
	return(optimizer)
}

optimization_constraint_weight<-function(
		optimizer,
		constraintType=c("=", ">=", "<="),
		constraintValue,
		symbols){
	util_validate(as.list(environment()))
	
	constraintMertic="POSITION_WEIGHT"
	optimizer=optimization_constraint(optimizer,constraintType,constraintMertic,constraintValue,symbols=symbols)
	return(optimizer)
}


util_optimizationFunction<-function(optimizer){
	balans=NULL
	portfolio=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",optimizer@portfolio),optimization_info=NULL)
	symbols=portfolio_symbols(portfolio)	
	rebalancingTimes=portfolio_value(portfolio)[,1]
	set=portfolio_getSettings(portfolio)
	set$resultsSamplingInterval='last'
	portfolio_settings(portfolio,set)

	for(time in rebalancingTimes){
		result<-.jcall(portfolio@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;", method="setToTime", util_POSIXTimeToDate(time))
		util_checkErrors(result)
		optimizer_temp<-optimization_goal(portfolio,goal=optimizer@goal, direction=optimizer@direction)
		if(!is.null(optimizer@portfolioValue)){
			optimization_constraint_portfolioValue(optimizer_temp,optimizer@portfolioValue)
		}
		for(i in 1:length(optimizer@functions)){
			function_names=names(formals(optimizer@functions[[i]]));
			
			temp_list=list()
			for(name in function_names){
				temp_list=c(temp_list,get(name))
			}
			names(temp_list)=function_names;
			constraintValue=tail(do.call(optimizer@functions[[1]],temp_list),1)
			optimizer_temp=optimization_constraint(optimizer_temp,optimizer@constraintTypeFunctions[[i]],optimizer@constraintMerticFunctions[[i]],as.double(constraintValue),optimizer@constraintConfidenceIntervalFunctions[[i]],optimizer@constraintSymbolsFunctions[[i]])
		}
			if(!is.null(optimizer@constraintMerticSimple)){
				for(j in 1:length(optimizer@constraintMerticSimple)){
					optimizer_temp=optimization_constraint(optimizer_temp,optimizer@constraintTypeSimple[[j]],optimizer@constraintMerticSimple[[j]],optimizer@constraintValueSimple[[j]],optimizer@constraintConfidenceInterval[[j]],optimizer@constraintSymbols[[j]])
				}
			}
		result<-.jcall(optimizer_temp@java,returnSig="Lcom/portfolioeffect/quant/client/result/MethodResult;",method="getOptimizedPortfolio")
		portfolio_TempOptim<-new("portfolio", java=getResult(result),optimization_info=NULL)
		temp=NULL
		for(symbol in symbols){
			temp=c(temp,position_quantity(portfolio_TempOptim,symbol)[2])
		}
		names(temp)=symbols;
		balans=cbind(balans,temp)
		
	}
	portfolioResult=new("portfolio", java=.jnew("com.portfolioeffect.quant.client.portfolio.Portfolio",optimizer@portfolio),optimization_info=NULL)
#	FromTime=.jcall(portfolioResult@java,returnSig="S",method="getFromTime")
	for(i in 1:length(symbols)){
		portfolio_addPosition(portfolioResult,symbols[i],balans[i,],time=util_POSIXTimeToDate(rebalancingTimes))
	}
	return(portfolioResult)
	
}

optimization_info<-function(portfolio){
	if(!is.null(portfolio@optimization_info)){
	nFunctionExecute=portfolio@optimization_info[1]
	nGlobalStart=portfolio@optimization_info[2]
	nLocalSolution=portfolio@optimization_info[3]
	nOptimizations=portfolio@optimization_info[4]
	nConstraintSatisfied=portfolio@optimization_info[5]
	
	r<-NULL
	r<-rbind(r,paste('Total number of calls to optimization goal method ',  	as.double(nFunctionExecute)))
	r<-rbind(r,paste('Total number of search paths of the optimization algorithm ',				as.double(nGlobalStart)))
	r<-rbind(r,paste('Total number of times all optimization constraints were satisfied ',					as.double(nConstraintSatisfied)))
	r<-rbind(r,paste('Total number of local optima found by the optimization algorithm ',		as.double(nLocalSolution)))
	r<-rbind(r,paste('Total number of optimizations ',			as.double(nOptimizations)))
	r<-rbind(r,paste('Average number of calls to optimization goal method per optimization step ',			round(as.double(nFunctionExecute)/as.double(nOptimizations),digits =3)))
	r<-rbind(r,paste('Average number of search paths of the optimization algorithm per optimization step ',	round(as.double(nGlobalStart)/as.double(nOptimizations),digits =3)))
	r<-rbind(r,paste('Average number of times all optimization constraints were satisfied per optimization step ',		round(as.double(nConstraintSatisfied)/as.double(nOptimizations),digits =3)))
	r<-rbind(r,paste('Average number of local optima found by the optimization algorithm per optimization step ',       		round(as.double(nLocalSolution)/as.double(nOptimizations),digits =3)))
	r<-rbind(r,paste('Total number of local optima skipped by the optimization algorithm ',        	as.double(nLocalSolution)-as.double(nOptimizations)))
	r<-rbind(r,paste('Average number of search path per local optimum (the lower, the better) ',	round(as.double(nGlobalStart)/as.double(nLocalSolution),digits =3)))
	r<-rbind(r,paste('Average number of local optima skipped by the optimization algorithm per optimization step ',	round((as.double(nLocalSolution)-as.double(nOptimizations))/as.double(nOptimizations),digits =3)))		
	
	colnames(r)<-array("",dim=1)
	rownames(r)<-array("",dim=NROW(r))
	cat (paste("OPTIMIZATION INFO",sep="    "))
	print(r, quote=FALSE)
}else{
	stop("given portfolio is not produced by optimization_run() method")
}
}