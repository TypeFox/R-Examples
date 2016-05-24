
########################################################################
#
#				Distributions
#
########################################################################

.makeDistributionDialog <- function(type="quantile"){

	#make distribution parameters
	meanParam <- new(ParamNumeric,"mean",0)

	sdParam <- new(ParamNumeric,"sd",1)
	sdParam$setLowerBound(0)

	meanlogParam <- new(ParamNumeric,"meanlog",0)

	sdlogParam <- new(ParamNumeric,"sdlog",1)
	sdlogParam$setLowerBound(0)

	dfParam <- new(ParamNumeric,"df")
	dfParam$setValue(1)
	dfParam$setLowerBound(0)

	df1Param <- new(ParamNumeric,"df1")
	df1Param$setValue(1)
	df1Param$setLowerBound(0)

	df2Param <- new(ParamNumeric,"df2")
	df2Param$setValue(1)
	df2Param$setLowerBound(0)

	ncpParam <- new(ParamNumeric,"ncp")
	ncpParam$setRequired(FALSE)

	rateParam <- new(ParamNumeric,"rate",1)
	rateParam$setLowerBound(0)

	maxParam <- new(ParamNumeric,"max",1)
	minParam <- new(ParamNumeric,"min",0)

	shapeParam <- new(ParamNumeric,"shape")

	shape1Param <- new(ParamNumeric,"shape1")
	shape1Param$setTitle("alpha")
	shape2Param <- new(ParamNumeric,"shape2")
	shape2Param$setTitle("beta")

	scaleParam <- new(ParamNumeric,"scale",1)
	locationParam <- new(ParamNumeric,"location",0)

	sizeParam <- new(ParamNumeric,"size")
	sizeParam$setLowerBound(0)
	probParam <- new(ParamNumeric,"prob")
	probParam$setLowerBound(0)
	probParam$setUpperBound(1)

	lambdaParam <- new(ParamNumeric,"lambda")
	lambdaParam$setLowerBound(0)

	mParam <- new(ParamNumeric,"m")
	mParam$setLowerBound(0)

	kParam <- new(ParamNumeric,"k")
	kParam$setLowerBound(0)
	
	xParam <- new(ParamNumeric,"x")
	qParam <- new(ParamNumeric,"q")
	pParam <- new(ParamNumeric,"p")
	pParam$setTitle("Probability")
	pParam$setLowerBound(0)
	pParam$setUpperBound(1)
	nParam <- new(ParamNumeric,"n")
	nParam$setLowerBound(0)

	lowerTailParam <- new(ParamLogical,"lower.tail",TRUE)

	#Define distributions
	addParams <-function(rFunc,params){
		for(i in 1:length(params))
			rFunc$add(params[[i]])
	}

	dnormFunc <- new(RFunction,"dnorm")
	addParams(dnormFunc,list(xParam,meanParam$clone(),sdParam$clone()))
	pnormFunc <- new(RFunction,"pnorm")
	addParams(pnormFunc,list(qParam,meanParam$clone(),sdParam$clone(),lowerTailParam$clone()))
	qnormFunc <- new(RFunction,"qnorm")
	addParams(qnormFunc,list(pParam,meanParam$clone(),sdParam$clone(),lowerTailParam$clone()))


	dtFunc <- new(RFunction,"dt")
	addParams(dtFunc,list(xParam,dfParam$clone(),ncpParam$clone()))
	ptFunc <- new(RFunction,"pt")
	addParams(ptFunc,list(qParam,dfParam$clone(),ncpParam$clone(),lowerTailParam$clone()))
	qtFunc <- new(RFunction,"qt")
	addParams(qtFunc,list(pParam,dfParam$clone(),ncpParam$clone(),lowerTailParam$clone()))

	ncp1Param <- new(ParamNumeric,"ncp",0)

	dchisqFunc <- new(RFunction,"dchisq")
	addParams(dchisqFunc,list(xParam,dfParam$clone(),ncp1Param$clone()))
	pchisqFunc <- new(RFunction,"pchisq")
	addParams(pchisqFunc,list(qParam,dfParam$clone(),ncp1Param$clone(),lowerTailParam$clone()))
	qchisqFunc <- new(RFunction,"qchisq")
	addParams(qchisqFunc,list(pParam,dfParam$clone(),ncp1Param$clone(),lowerTailParam$clone()))


	dfFunc <- new(RFunction,"df")
	addParams(dfFunc,list(xParam,df1Param$clone(),df2Param$clone(),ncpParam$clone()))
	pfFunc <- new(RFunction,"pf")
	addParams(pfFunc,list(qParam,df1Param$clone(),df2Param$clone(),ncpParam$clone(),lowerTailParam$clone()))
	qfFunc <- new(RFunction,"qf")
	addParams(qfFunc,list(pParam,df1Param$clone(),df2Param$clone(),ncpParam$clone(),lowerTailParam$clone()))


	dexpFunc <- new(RFunction,"dexp")
	addParams(dexpFunc,list(xParam,rateParam$clone()))
	pexpFunc <- new(RFunction,"pexp")
	addParams(pexpFunc,list(qParam,rateParam$clone(),lowerTailParam$clone()))
	qexpFunc <- new(RFunction,"qexp")
	addParams(qexpFunc,list(pParam,rateParam$clone(),lowerTailParam$clone()))


	dunifFunc <- new(RFunction,"dunif")
	addParams(dunifFunc,list(xParam,minParam$clone(),maxParam$clone()))
	punifFunc <- new(RFunction,"punif")
	addParams(punifFunc,list(qParam,minParam$clone(),maxParam$clone(),lowerTailParam$clone()))
	qunifFunc <- new(RFunction,"qunif")
	addParams(qunifFunc,list(pParam,minParam$clone(),maxParam$clone(),lowerTailParam$clone()))


	dbetaFunc <- new(RFunction,"dbeta")
	addParams(dbetaFunc,list(xParam,shape1Param$clone(),shape2Param$clone()))
	pbetaFunc <- new(RFunction,"pbeta")
	addParams(pbetaFunc,list(qParam,shape1Param$clone(),shape2Param$clone(),lowerTailParam$clone()))
	qbetaFunc <- new(RFunction,"qbeta")
	addParams(qbetaFunc,list(pParam,shape1Param$clone(),shape2Param$clone(),lowerTailParam$clone()))


	dcauchyFunc <- new(RFunction,"dcauchy")
	addParams(dcauchyFunc,list(xParam,locationParam$clone(),scaleParam$clone()))
	pcauchyFunc <- new(RFunction,"pcauchy")
	addParams(pcauchyFunc,list(qParam,locationParam$clone(),scaleParam$clone(),lowerTailParam$clone()))
	qcauchyFunc <- new(RFunction,"qcauchy")
	addParams(qcauchyFunc,list(pParam,locationParam$clone(),scaleParam$clone(),lowerTailParam$clone()))


	dlogisFunc <- new(RFunction,"dlogis")
	addParams(dlogisFunc,list(xParam,locationParam$clone(),scaleParam$clone()))
	plogisFunc <- new(RFunction,"plogis")
	addParams(plogisFunc,list(qParam,locationParam$clone(),scaleParam$clone(),lowerTailParam$clone()))
	qlogisFunc <- new(RFunction,"qlogis")
	addParams(qlogisFunc,list(pParam,locationParam$clone(),scaleParam$clone(),lowerTailParam$clone()))


	dlnormFunc <- new(RFunction,"dlnorm")
	addParams(dlnormFunc,list(xParam,meanlogParam$clone(),sdlogParam$clone()))
	plnormFunc <- new(RFunction,"plnorm")
	addParams(plnormFunc,list(qParam,meanlogParam$clone(),sdlogParam$clone(),lowerTailParam$clone()))
	qlnormFunc <- new(RFunction,"qlnorm")
	addParams(qlnormFunc,list(pParam,meanlogParam$clone(),sdlogParam$clone(),lowerTailParam$clone()))


	dgammaFunc <- new(RFunction,"dgamma")
	addParams(dgammaFunc,list(xParam,shapeParam$clone(),rateParam$clone()))
	pgammaFunc <- new(RFunction,"pgamma")
	addParams(pgammaFunc,list(qParam,shapeParam$clone(),rateParam$clone(),lowerTailParam$clone()))
	qgammaFunc <- new(RFunction,"qgamma")
	addParams(qgammaFunc,list(pParam,shapeParam$clone(),rateParam$clone(),lowerTailParam$clone()))


	dweibullFunc <- new(RFunction,"dweibull")
	addParams(dweibullFunc,list(xParam,shapeParam$clone(),scaleParam$clone()))
	pweibullFunc <- new(RFunction,"pweibull")
	addParams(pweibullFunc,list(qParam,shapeParam$clone(),scaleParam$clone(),lowerTailParam$clone()))
	qweibullFunc <- new(RFunction,"qweibull")
	addParams(qweibullFunc,list(pParam,shapeParam$clone(),scaleParam$clone(),lowerTailParam$clone()))


	dbinomFunc <- new(RFunction,"dbinom")
	addParams(dbinomFunc,list(xParam,sizeParam$clone(),probParam$clone()))
	pbinomFunc <- new(RFunction,"pbinom")
	addParams(pbinomFunc,list(qParam,sizeParam$clone(),probParam$clone(),lowerTailParam$clone()))
	qbinomFunc <- new(RFunction,"qbinom")
	addParams(qbinomFunc,list(pParam,sizeParam$clone(),probParam$clone(),lowerTailParam$clone()))


	dpoisFunc <- new(RFunction,"dpois")
	addParams(dpoisFunc,list(xParam,lambdaParam$clone()))
	ppoisFunc <- new(RFunction,"ppois")
	addParams(ppoisFunc,list(qParam,lambdaParam$clone(),lowerTailParam$clone()))
	qpoisFunc <- new(RFunction,"qpois")
	addParams(qpoisFunc,list(pParam,lambdaParam$clone(),lowerTailParam$clone()))


	dgeomFunc <- new(RFunction,"dgeom")
	addParams(dgeomFunc,list(xParam,probParam$clone()))
	pgeomFunc <- new(RFunction,"pgeom")
	addParams(pgeomFunc,list(qParam,probParam$clone(),lowerTailParam$clone()))
	qgeomFunc <- new(RFunction,"qgeom")
	addParams(qgeomFunc,list(pParam,probParam$clone(),lowerTailParam$clone()))


	dhyperFunc <- new(RFunction,"dhyper")
	addParams(dhyperFunc,list(xParam,mParam$clone(),nParam$clone(),kParam$clone()))
	phyperFunc <- new(RFunction,"phyper")
	addParams(phyperFunc,list(qParam,mParam$clone(),nParam$clone(),kParam$clone(),lowerTailParam$clone()))
	qhyperFunc <- new(RFunction,"qhyper")
	addParams(qhyperFunc,list(pParam,mParam$clone(),nParam$clone(),kParam$clone(),lowerTailParam$clone()))

	
	panelView <- new(ParamNumeric,"")$VIEW_RFUNCTION_PANEL
	
	#make function list and add functions
	quantileFuncList <- new(RFunctionList,"Quantile function");

	#parameters common to all functions should go at the top
	globals <- .jarray(list(pParam),"org.rosuda.deducer.widgets.param.Param")
	quantileFuncList$setGlobalParams(globals)

	quantileFuncList$setViewType(panelView)
	quantileFuncList$addRFunction("Normal",qnormFunc)
	quantileFuncList$addRFunction("t",qtFunc)
	quantileFuncList$addRFunction("Chi-squared",qchisqFunc)
	quantileFuncList$addRFunction("F",qfFunc)
	quantileFuncList$addRFunction("Exponential",qexpFunc)
	quantileFuncList$addRFunction("Uniform",qunifFunc)
	quantileFuncList$addRFunction("Beta",qbetaFunc)
	quantileFuncList$addRFunction("Cauchy",qcauchyFunc)
	quantileFuncList$addRFunction("Logistic",qlogisFunc)
	quantileFuncList$addRFunction("Log normal",qlnormFunc)
	quantileFuncList$addRFunction("Gamma",qgammaFunc)
	quantileFuncList$addRFunction("Weibull",qweibullFunc)
	quantileFuncList$addRFunction("Binomial",qbinomFunc)
	quantileFuncList$addRFunction("Poisson",qpoisFunc)
	quantileFuncList$addRFunction("Geometric",qgeomFunc)
	quantileFuncList$addRFunction("Hyper geometric",qhyperFunc)


	#make function list and add functions
	distFuncList <- new(RFunctionList,"Distribution function");

	#parameters common to all functions should go at the top
	globals <- .jarray(list(xParam),"org.rosuda.deducer.widgets.param.Param")
	distFuncList$setGlobalParams(globals)

	distFuncList$setViewType(panelView)
	distFuncList$addRFunction("Normal",dnormFunc)
	distFuncList$addRFunction("t",dtFunc)
	distFuncList$addRFunction("Chi-squared",dchisqFunc)
	distFuncList$addRFunction("F",dfFunc)
	distFuncList$addRFunction("Exponential",dexpFunc)
	distFuncList$addRFunction("Uniform",dunifFunc)
	distFuncList$addRFunction("Beta",dbetaFunc)
	distFuncList$addRFunction("Cauchy",dcauchyFunc)
	distFuncList$addRFunction("Logistic",dlogisFunc)
	distFuncList$addRFunction("Log normal",dlnormFunc)
	distFuncList$addRFunction("Gamma",dgammaFunc)
	distFuncList$addRFunction("Weibull",dweibullFunc)
	distFuncList$addRFunction("Binomial",dbinomFunc)
	distFuncList$addRFunction("Poisson",dpoisFunc)
	distFuncList$addRFunction("Geometric",dgeomFunc)
	distFuncList$addRFunction("Hyper geometric",dhyperFunc)


	#make function list and add functions
	cumDistFuncList <- new(RFunctionList,"Cumulative distribution function");

	#parameters common to all functions should go at the top
	globals <- .jarray(list(qParam),"org.rosuda.deducer.widgets.param.Param")
	cumDistFuncList$setGlobalParams(globals)

	cumDistFuncList$setViewType(panelView)
	cumDistFuncList$addRFunction("Normal",pnormFunc)
	cumDistFuncList$addRFunction("t",ptFunc)
	cumDistFuncList$addRFunction("Chi-squared",pchisqFunc)
	cumDistFuncList$addRFunction("F",pfFunc)
	cumDistFuncList$addRFunction("Exponential",pexpFunc)
	cumDistFuncList$addRFunction("Uniform",punifFunc)
	cumDistFuncList$addRFunction("Beta",pbetaFunc)
	cumDistFuncList$addRFunction("Cauchy",pcauchyFunc)
	cumDistFuncList$addRFunction("Logistic",plogisFunc)
	cumDistFuncList$addRFunction("Log normal",plnormFunc)
	cumDistFuncList$addRFunction("Gamma",pgammaFunc)
	cumDistFuncList$addRFunction("Weibull",pweibullFunc)
	cumDistFuncList$addRFunction("Binomial",pbinomFunc)
	cumDistFuncList$addRFunction("Poisson",ppoisFunc)
	cumDistFuncList$addRFunction("Geometric",pgeomFunc)
	cumDistFuncList$addRFunction("Hyper geometric",phyperFunc)

	if(type=="quantile")
		model <- quantileFuncList
	else if(type=="distribution")
		model <- distFuncList
	else if(type=="CDF")
		model <- cumDistFuncList
	
	rfd <- new(RFunctionListDialog, model )
	rfd$setSize(475L,700L)
	rfd$setLocationRelativeTo(.jnull())
	
	rfd
}







