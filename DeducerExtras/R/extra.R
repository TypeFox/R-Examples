# TODO: Add comment
# 
# Author: Ian
###############################################################################


########################################################################
#
#				Summarize data.frame
#
########################################################################
.makeSummaryDialog <- function(){
	summaryFunc <- new(RFunction,"summary")
	
	dataParam <- new(ParamRObject)
	dataParam$setRObjectClass("data.frame")
	dataParam$setName(.jnull())
	dataParam$setTitle("data")
	summaryFunc$add(dataParam)
	
	rfd <- new(RFunctionDialog, summaryFunc)
	rfd$setSize(250L,150L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}

########################################################################
#
#				Load data from package
#
########################################################################
.makePackageDataDialog <- function(){
	dat <- data()$results
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamCharacter <- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	rf <- new(RFunction,"data")
	rf$setTitle("Load Data From Package")
	
	p <- new(ParamCharacter,"data")
	p$setTitle("Data set")
	p$setOptions(dat[,3])
	p$setLabels(dat[,4])
	p$setViewType(p$VIEW_COMBO)
	rf$add(p)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setSize(250L,125L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}



########################################################################
#
#				Single proportion: asymptotic
#
########################################################################

.makeProportionDialog <- function(){
	oneSample <- new(RFunction,"one.sample.test")
	oneSample$setRequiresVariableSelector(TRUE)
	
	datParam <- new(ParamMultipleVariables,"variables")
	datParam$setFormat(datParam$FORMAT_VARIABLE)
	datParam$setTitle("Variables")
	oneSample$add(datParam)
	
	p <- new(ParamNumeric,"p")
	p$setValue(.5)
	p$setLowerBound(0)
	p$setUpperBound(1)
	oneSample$add(p)
	
	alt <- new(ParamCharacter,"alternative","two.sided")
	alt$setTitle("Type")
	alt$setOptions(c("two.sided","less","greater"))
	alt$setViewType(alt$VIEW_COMBO)
	oneSample$add(alt)
	
	
	func <- "func<-function(x,...) {
			y<-as.factor(na.omit(x))
			if(length(levels(y))>2) stop()
			prop.test(sum(y == rev(levels(y))[1]),n=length(y),...)
			}"
	test <- new(ParamAny,"test")
	test$setViewType(test$VIEW_HIDDEN)
	test$setValue(func)
	oneSample$add(test)
	
	conf <- new(ParamNumeric,"conf.level",.95)
	conf$setLowerBound(0)
	conf$setUpperBound(1)
	oneSample$add(conf)
	
	corr <- new(ParamLogical,"correct",FALSE)
	corr$setDefaultValue(TRUE)
	corr$setViewType(corr$VIEW_HIDDEN)
	oneSample$add(corr)
	
	rfd <- new(RFunctionDialog, oneSample)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


########################################################################
#
#				Single proportion: exact
#
########################################################################

.makeExactProportionDialog <- function(){
	oneSample <- new(RFunction,"one.sample.test")
	oneSample$setRequiresVariableSelector(TRUE)
	
	datParam <- new(ParamMultipleVariables,"variables")
	datParam$setFormat(datParam$FORMAT_VARIABLE)
	datParam$setTitle("Variables")
	oneSample$add(datParam)
	
	p <- new(ParamNumeric,"p")
	p$setValue(.5)
	p$setLowerBound(0)
	p$setUpperBound(1)
	oneSample$add(p)
	
	alt <- new(ParamCharacter,"alternative","two.sided")
	alt$setTitle("Type")
	alt$setOptions(c("two.sided","less","greater"))
	alt$setViewType(alt$VIEW_COMBO)
	oneSample$add(alt)
	
	
	func <- "func<-function(x,...) {
			y<-as.factor(na.omit(x))
			if(length(levels(y))>2) stop()
			binom.test(sum(y == rev(levels(y))[1]),n=length(y),...)
			}"
	test <- new(ParamAny,"test")
	test$setViewType(test$VIEW_HIDDEN)
	test$setValue(func)
	oneSample$add(test)
	
	conf <- new(ParamNumeric,"conf.level",.95)
	conf$setLowerBound(0)
	conf$setUpperBound(1)
	oneSample$add(conf)
	
	rfd <- new(RFunctionDialog, oneSample)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}



########################################################################
#
#				n proportions: asymptotic
#
########################################################################
.makeNProportionDialog <- function(){
	nSample <- new(RFunction,"prop.test")
	nSample$setRequiresVariableSelector(TRUE)
	
	tableFunc <- new(RFunction,"table")
	tableFunc$setRequiresVariableSelector(TRUE)
	nSample$add(tableFunc)
	
	var1 <- new(ParamVariable,"variable")
	var2 <- new(ParamVariable,"group")
	tableFunc$add(var2)
	tableFunc$add(var1)
	
	alt <- new(ParamCharacter,"alternative","two.sided")
	alt$setTitle("Type")
	alt$setOptions(c("two.sided","less","greater"))
	alt$setViewType(alt$VIEW_COMBO)
	nSample$add(alt)
	
	conf <- new(ParamNumeric,"conf.level",.95)
	conf$setLowerBound(0)
	conf$setUpperBound(1)
	nSample$add(conf)
	
	corr <- new(ParamLogical,"correct",FALSE)
	corr$setDefaultValue(TRUE)
	corr$setViewType(corr$VIEW_HIDDEN)
	nSample$add(corr)
	
	rfd <- new(RFunctionDialog, nSample)
	rfd$setSize(600L,400L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


########################################################################
#
#				k-sample Equal variance
#
########################################################################
.makeEqualVarianceDialog <- function(){
	levene <- new(RFunction,"levene.test")
	bartlett <- new(RFunction, "bartlett.test")
	
	var1 <- new(ParamVariable)
	var1$setTitle("Variable")
	var2 <- new(ParamVariable)
	var2$setTitle("Group")
	levene$add(var1)
	levene$add(var2)
	bartlett$add(var1)
	bartlett$add(var2)
	
	opt <- new(ParamCharacter,"option","mean")
	opt$setTitle("Method")
	opt$setOptions(c("mean", "median","trim.mean"))
	opt$setLabels(c("Levene (mean)", "Brown-Forsythe (median)","Trimmed mean"))
	opt$setViewType(opt$VIEW_COMBO)
	levene$add(opt)
	
	trim <- new(ParamNumeric,"trim.alpha")
	trim$setTitle("% trimmed from each tail")
	trim$setRequired(FALSE)
	trim$setLowerBound(0)
	trim$setUpperBound(.5)
	levene$add(trim)
	
	#make function list and add functions
	prf <- new(RFunctionList,"K-sample equality of variance");
	prf$setViewType(prf$VIEW_RFUNCTION_PANEL)
	prf$addRFunction("Levene's test",levene,FALSE)
	prf$addRFunction("Bartlett's test",bartlett,FALSE)
	prf$setRequiresVariableSelector(TRUE)
	prf$setActiveFunctions("Levene's test")
	#parameters common to all functions should go at the top
	globals <- .jarray(list(var1,var2),"org.rosuda.deducer.widgets.param.Param")
	prf$setGlobalParams(globals)
	
	
	#make dialog and display
	rfd <- new(RFunctionListDialog, prf )
	rfd$setSize(500L,570L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}

########################################################################
#
#				t-test power
#
########################################################################
.makeTTestPowerDialog <- function(){
	status <- .deducer$requirePackage("pwr")
	if(status=="installed"){
		execute("library('pwr')")
	}else if(status=="not-installed"){
		stop("package pwr required")
	}
	dialog <- new(SimpleRDialog)
	dialog$setSize(350L,580L)
	dialog$setTitle("t-test power analysis")
	
	#type of test
	test<- new(ComboBoxWidget,"type of test",c("two.sample", "one.sample", "paired"))
	test$setDefaultModel("two.sample")
	addComponent(dialog, test,100,900,200, 100)
	
	#Sample size
	ss <- new(TextAreaWidget,"Sample size")
	addComponent(dialog, ss,210,700,310, 300)
	
	#sig
	sig <- new(TextAreaWidget,"significance level")
	sig$setDefaultModel("0.05")
	addComponent(dialog, sig,320,700,420, 300)
	
	#power
	pow <- new(TextAreaWidget,"Power")
	pow$setDefaultModel("0.80")
	addComponent(dialog, pow,430,700,530, 300)
	
	#effect size
	eff <- new(TextAreaWidget,"Cohens D")
	eff$setDefaultModel(".5")
	addComponent(dialog, eff,540,700,640, 300)
	
	
	#alternative
	test<- new(ComboBoxWidget,"alternative",c("two.sided", "less","greater"))
	test$setDefaultModel("two.sided")
	addComponent(dialog, test,650,900,750, 100)
	
	
	runDialog <- function(state){
		#print(state)
		cmd <- "require(pwr)\npwr.t.test("
		
		if(state[['Sample size']] == "")
			parameter <- "n=NULL"
		else
			parameter = paste("n=",state[['Sample size']],sep="")
		cmd <- paste(cmd,parameter);
		
		if(state[['significance level']] == "")
			parameter <- ",sig.level=NULL"
		else
			parameter = paste(",sig.level=",state[['significance level']],sep="")
		cmd <- paste(cmd,parameter);
		
		if(state[['Power']] == "")
			parameter <- ",power=NULL"
		else
			parameter = paste(",power=",state[['Power']],sep="")
		cmd <- paste(cmd,parameter);
		
		if(state[['Cohens D']] == "")
			parameter <- ",d=NULL"
		else
			parameter = paste(",d=",state[['Cohens D']],sep="")
		cmd <- paste(cmd,parameter);
		
		parameter = paste(",alternative='",state[['alternative']],"'",sep="")
		cmd <- paste(cmd,parameter);
		
		parameter = paste(",type='",state[['type of test']],"')",sep="")
		cmd <- paste(cmd,parameter);
		
		execute(cmd)
	}
	
	dialog$setRunFunction(toJava(runDialog))
	dialog
}
########################################################################
#
#				k-means cluster
#
########################################################################
.makeKMeansDialog <- function(){
	naOmitFunc <- new(RFunction,"na.omit")
	naOmitFunc$setRequiresVariableSelector(TRUE)
	
	datParam <- new(ParamMultipleVariables,"object")
	datParam$setTitle("Variables")
	naOmitFunc$add(datParam)
	
	kmeansFunc <- new(RFunction,"kmeans")
	kmeansFunc$setTitle("Clustering")
	kmeansFunc$add(naOmitFunc)
	
	centersParam <- new(ParamNumeric,"centers")
	centersParam$setTitle("# of clusters")
	centersParam$setValue(2)
	centersParam$setLowerBound(1)
	kmeansFunc$add(centersParam)
	
	
	algParam <- new(ParamCharacter,"algorithm","Hartigan-Wong")
	algParam$setTitle("Type")
	algParam$setOptions(c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen"))
	algParam$setViewType(algParam$VIEW_COMBO)
	kmeansFunc$add(algParam)
	
	centersParam <- new(ParamNumeric,"iter.max",10)
	centersParam$setTitle("maximum # of iterations")
	centersParam$setLowerBound(1)
	kmeansFunc$add(centersParam)
	
	clustFuncList <- new(RFunctionList,"k-means clustring");
	clustFuncList$setViewType(clustFuncList$VIEW_RFUNCTION_PANEL)
	clustFuncList$setRequiresVariableSelector(TRUE)
	clustFuncList$addRFunction("Cluster",kmeansFunc,TRUE,TRUE,TRUE,"kmeansResult")
	
	
	#make dialog and display
	rfd <- new(RFunctionListDialog, clustFuncList )
	rfd$setSize(550L,500L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}



########################################################################
#
#				Apply k-means
#
########################################################################

predict.kmeans <- function(object,data=NULL,...){
	if(is.null(data))
		return(object$cluster)
	centers <- object$centers
	vars <- colnames(centers)
	dat <- data[,vars,drop=FALSE]
	clusters <- rep(NA,nrow(dat))
	for(i in 1:nrow(dat)){
		obs <- dat[i,]
		dists <- apply(centers,1,function(x) dist(rbind(obs,x)))
		clust <- names(which.min(dists))
		if(length(clust)>0)
			clusters[i] <- clust
		else
			clusters[i] <- NA
		
	}
	as.numeric(clusters)
}

applyModel <- function(object,data,...) data.frame(data,predict(object,data=data,...))

.makeApplyKMeansDialog <- function(){

	assignFunc <- new(RFunction,"assign")
	assignFunc$setTitle("apply k-means")
	
	newDataParam <- new(ParamCharacter)
	newDataParam$setTitle("generated data name:")
	newDataParam$setValue("kmeansData")
	newDataParam$setName(.jnull())
	
	applyFunc <- new(RFunction,"applyModel")
	
	modelParam <- new(ParamRObject,"object")
	modelParam$setTitle("model")
	modelParam$setRObjectClass("kmeans")
	
	dataParam <- new(ParamRObject,"data")
	dataParam$setRObjectClass("data.frame")
	
	assignFunc$add(newDataParam)
	assignFunc$add(applyFunc)
	applyFunc$add(dataParam)
	applyFunc$add(modelParam)
	
	rfd <- new(RFunctionDialog,assignFunc)
	rfd$setSize(450L,250L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


########################################################################
#
#				hierarchical cluster
#
########################################################################

.makeHClustDialog <- function(){

	clustFuncList <- new(RFunctionList,"Hierarchical clustring");
	clustFuncList$setViewType(clustFuncList$VIEW_RFUNCTION_PANEL)
	clustFuncList$setRequiresVariableSelector(TRUE)
	
	distFunc <- new(RFunction,"dist")
	
	x <- new(ParamMultipleVariables,"x")
	x$setTitle("Data")
	
	method <- new(ParamCharacter,"method","euclidean")
	method$setTitle("Distance")
	method$setOptions(c("euclidean", "maximum","manhattan","canberra","binary","minkowski"))
	method$setViewType(method$VIEW_COMBO)
	
	hclustFunc <- new(RFunction,"hclust")
	
	hcMethod <- new(ParamCharacter,"method","complete")
	hcMethod$setTitle("Method")
	hcMethod$setOptions(c("ward", "single", "complete", 
					"average", "mcquitty", "median" , "centroid"))
	hcMethod$setViewType(hcMethod$VIEW_COMBO)
	
	dd <- new(ParamRFunctionResult,clustFuncList,"dist")
	dd$setName("d")
	
	plotFunc <- new(RFunction,"plot")
	
	clustResult <- new(ParamRFunctionResult,clustFuncList,"hclust")
	clustResult$setName(.jnull())
	
	cutreeFunc <- new(RFunction,"cutree")
	cutreeFunc$setTitle("Cluster membership")
	
	k <- new(ParamNumeric,"k")
	k$setTitle("# of clusters")
	k$setLowerBound(0)
	
	distFunc$add(x)
	distFunc$add(method)
	hclustFunc$add(dd)
	hclustFunc$add(hcMethod)
	plotFunc$add(clustResult)
	cutreeFunc$add(clustResult)
	cutreeFunc$add(k)
	
	clustFuncList$addRFunction("dist",distFunc,TRUE,FALSE,FALSE,"<auto>")
	clustFuncList$addRFunction("hclust",hclustFunc,TRUE,TRUE,TRUE,"<auto>")
	clustFuncList$addRFunction("Dendogram",plotFunc,FALSE,FALSE,FALSE,"<auto>")
	clustFuncList$addRFunction("Cluster groups",cutreeFunc)
	#make dialog and display
	rfd <- new(RFunctionListDialog, clustFuncList )
	rfd$setSize(470L,620L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}



########################################################################
#
#				Multi dimensional scaling
#
########################################################################


.makeMDSDialog <- function(){

	mdsFuncList <- new(RFunctionList,"Multi-dimensional scaling");
	mdsFuncList$setViewType(mdsFuncList$VIEW_RFUNCTION_PANEL)
	mdsFuncList$setRequiresVariableSelector(TRUE)
	
	distFunc <- new(RFunction,"dist")
	
	x <- new(ParamMultipleVariables,"x")
	x$setTitle("Data")
	
	method <- new(ParamCharacter,"method","euclidean")
	method$setTitle("Distance")
	method$setOptions(c("euclidean", "maximum","manhattan","canberra","binary","minkowski"))
	method$setViewType(method$VIEW_COMBO)
	
	dd <- new(ParamRFunctionResult,mdsFuncList,"Distance")
	dd$setName(.jnull())
	
	mdsFunc <- new(RFunction,"cmdscale")
	mdsFunc$setTitle("Scaling")
	
	res <- new(ParamRFunctionResult,mdsFuncList,"Scaling")
	res$setName(.jnull())
	
	k <- new(ParamNumeric,"k")
	k$setTitle("# of dimensions")
	k$setValue(2)
	k$setLowerBound(0)
	
	plotFunc <- new(RFunction,"plot")
	
	distFunc$add(x)
	distFunc$add(method)
	mdsFunc$add(dd)
	mdsFunc$add(k)
	plotFunc$add(res)
	
	mdsFuncList$addRFunction("Distance",distFunc,TRUE,FALSE,FALSE,"<auto>")
	mdsFuncList$addRFunction("Scaling",mdsFunc,TRUE,FALSE,TRUE,"<auto>")
	mdsFuncList$addRFunction("Plot first two dimensions",plotFunc,FALSE,FALSE,FALSE,"<auto>")
	#make dialog and display
	rfd <- new(RFunctionListDialog, mdsFuncList )
	rfd$setSize(520L,550L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}




########################################################################
#
#				Ranking analysis
#
########################################################################


.makeRankingAnalysisDialog <- function(){
	# dialog frame
	dialog <- new(SimpleRDialog)
	dialog$setSize(400L,400L)
	dialog$setTitle("Ranking analysis")
	# add list of current variables
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(dialog,variableSelector,25,400,850,25)
	# add list of selected variables
	variableList <- new(VariableListWidget,variableSelector)
	variableList$setTitle("variables")
	addComponent(dialog,variableList,100,950,550,420)
	variableListLabel <- new(JLabel,"selected")
	addComponent(dialog,variableListLabel,40,950,80,600)
	# add analysis options
	testBoxes <- new(CheckBoxesWidget,"Tests",c("Friedman","Kendall-W"))
	addComponent(dialog,testBoxes,600,950,850,550)
	# add help button
	dialog$addHelpButton("pmwiki.php?n=Main.RankingAnalysis")

	# Analysis Check function
	.rankingAnalysisCheck <- function(state){
		# ensure that at least two variables are selected
		if (length(state$variables)<2)
			return("At least two variables must be selected")
		return("")
	}
	
	# Run function
	.rankingAnalysisRun <- function(state){
		# create object for ranking analysis
		cmd <- "ranking.analysis <- list()\n "
		# make ranking matrix from data frame
		# first copy the values from the data frame
		cmd <- paste(cmd,"ranking.analysis$rmat <- cbind(",state$variables[1],"=",state$data,"$",state$variables[1],sep="")
		for (varname in state$variables[-1]){
			cmd <- paste(cmd,",",varname,"=",state$data,"$",varname,sep="")
			}
		cmd <- paste(cmd,")\n")
		# then convert the rows into rankings
		cmd <- paste(cmd,"for (judge in 1:nrow(ranking.analysis$rmat)) ranking.analysis$rmat[judge,] <- rank(ranking.analysis$rmat[judge,])\n")
		cmd <- paste(cmd,"rm(judge)\n")
		# print summary of rankings
		cmd <- paste(cmd,"summary(ranking.analysis$rmat)\n")
		# Friedman test (if selected)
		if ("Friedman" %in% state$Tests){
			cmd <- paste(cmd,"ranking.analysis$friedman.results <- friedman.test(ranking.analysis$rmat)\n")
			cmd <- paste(cmd,"ranking.analysis$friedman.results\n")
		}
		# Kendall-W test (if selected)
		if ("Kendall-W" %in% state$Tests){
			require(irr)
			cmd <- paste(cmd,"ranking.analysis$kendall.results <- kendall(t(ranking.analysis$rmat),correct=TRUE)\n")	
			cmd <- paste(cmd,"ranking.analysis$kendall.results\n")
		}
		# execute command
		execute(cmd)
	}


	# check and run functions (defined below)
	dialog$setCheckFunction(toJava(.rankingAnalysisCheck))
	dialog$setRunFunction(toJava(.rankingAnalysisRun))
	# run dialog
	return(dialog)
}


########################################################################
#
#				3d scatter plot
#
########################################################################

.make3dScatterPlotDialog <- function(){
	fun <- new(RFunction,"scatter3d")
	
	if(!require("rgl")){
		install.packages("rgl")
		library("rgl")
	}
	
	var1 <- new(ParamVariable)
	var1$setTitle("x")
	var2 <- new(ParamVariable)
	var2$setTitle("y")
	var3 <- new(ParamVariable)
	var3$setTitle("z")
	var4 <- new(ParamVariable,"group")
	var4$setRequired(FALSE)
	lg <- new(ParamLogical)
	lg$setTitle("Show surface")
	lg$setName("surface")
	lg$setDefaultValue(TRUE)
	lg$setValue(FALSE)
	
	fun$add(var1)
	fun$add(var2)
	fun$add(var3)
	fun$add(var4)
	fun$add(lg)
	
	rfd <- new(RFunctionDialog, fun)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


########################################################################
#
#				Vignette browser
#
########################################################################

# Create list of vignettes
.makeVignetteDialog <- function(){
	vig.results <- vignette()$results[,c("Title","Package","Item")]
	vig.titles <- gsub("\\(source, pdf\\)|\\\"|\\\\","",vig.results[,1])
	vig.labels <- paste0(vig.titles," (",vig.results[,2],")")
	# Make dialog
	dialog <- new(SimpleRDialog)
	dialog$setTitle("Open Vignette")
	dialog$setSize(400L,150L)
	# Add vignette selector
	vignetteSelector <- new(ComboBoxWidget,"Vignettes",vig.labels)
	addComponent(dialog,vignetteSelector,50L,950L,500L,50L)
	# Run function
	vignetteRunFunction <- function(state){
		sel <- which(vig.labels == state)
		cmd <- paste0('vignette("', vig.results[sel[1],3], '", package="', vig.results[sel[1],2], '")')
		execute(cmd)
	}
	# Change OkayCancel buttons (remove Reset)
	dialog$setOkayCancel(FALSE,TRUE,dialog)
	dialog$setRunFunction(toJava(vignetteRunFunction))
	dialog
}



