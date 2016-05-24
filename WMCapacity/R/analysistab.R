
.womActiveAnalysisTab<-function(status=TRUE)
{

	theWidget("analysisPageBox")$setSensitive(status)

}

.womInitializeAnalysisPageColumns <- function()
{
	.womCreateAnalysisDefinedModelsColumns()
}

# adapted from RGtk2 treeStore demo
.womToggleAnalyzeModel<- function(cell, path.str, data)
{
  #checkPtrType(data, "GtkTreeModel")
  #model <- theWidget("analysisDefinedModelsTreeview")$getModel()
  model <- theWidget("modelsDefinedModelsTreeview")$getModel()

  #model2 <- gtkTreeModelSortGetModel(model)
 
  path <- gtkTreePathNewFromString(path.str)
  
  #column <- cell$getData("column")
  column <- .womDefinedModelsTreeCols("toBeAnalyzed")
  
  # get toggled iter
  iter <- model$getIter(path)$iter
  toggle.item <- gtkTreeModelGetValue(model, iter, column)$value
  
  # do something with the value
  toggle.item <- !toggle.item

  # set new value
  model$set(iter, column, toggle.item)
  #gtkTreeStoreSetValue(model2, iter, column, toggle.item)
}

.womCreateAnalysisDefinedModelsColumns <- function()
{
	#modelsTreeview <- theWidget("modelsDefinedModelsTreeview")
	analysisTreeview <- theWidget("analysisDefinedModelsTreeview")
	
	model = theWidget("modelsDefinedModelsTreeview")$getModel()
	#gtkTreeViewSetModel(analysisTreeview, model)	
	
	
	treeview = analysisTreeview
	#Create effects tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
	
		# To be analyzed
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Analyze", renderer, 
  								active = .womDefinedModelsTreeCols("toBeAnalyzed")
								)
		
		gSignalConnect(renderer, "toggled", .womToggleAnalyzeModel, model)
		
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

	
		# name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Name", renderer, 
  								text = .womDefinedModelsTreeCols("name"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# hasResults
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Results?", renderer, 
  								active = .womDefinedModelsTreeCols("hasResults")
								)
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(FALSE)
		
		
		# time
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Time", renderer, 
  								text = .womDefinedModelsTreeCols("timeAnalyzed"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

				
		# acc rate
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Acc. Rate", renderer, 
  								text = .womDefinedModelsTreeCols("accRate"),
								foreground = .womDefinedModelsTreeCols("accColor"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column, .womDefinedModelsTreeCols("accRate"))


		# iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("iterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# effective iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effective Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("effectiveIterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# burnin
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Burnin", renderer, 
  								text = .womDefinedModelsTreeCols("burnin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		# hybrid EPS 
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Lower eps.", renderer, 
  								text = .womDefinedModelsTreeCols("epsLow"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Upper eps.", renderer, 
  								text = .womDefinedModelsTreeCols("epsUpp"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# Leapfrog	
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Leapfrog", renderer, 
  								text = .womDefinedModelsTreeCols("leapfrog"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# use Metropolis
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH instead?", renderer, 
  								active = .womDefinedModelsTreeCols("useMet")
								)
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(FALSE)
	
		# metrop scale
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH scale", renderer, 
  								text = .womDefinedModelsTreeCols("metropScale"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# metrop scale
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH thin", renderer, 
  								text = .womDefinedModelsTreeCols("metropThin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}

}



.clicked_analysis_load_selected <- function(treeview, path, view_column)
{
	treeview <- theWidget("analysisDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetAnalysisModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	
	.womLoadModelMCMCSettings(modelName)
}

.womGetAnalysisModelSelection<-function()
{
	treeview = theWidget("analysisDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	#vis <- model$get(iter, .womDefinedModelsTreeCols("visible"))[[1]]
	return(list(iter=iter,visible=NULL))
}

.womLoadModelMCMCSettings <- function(modelName)
{
	myModel <- wommbatAnalysis$Models[[modelName]]

	if(is.null(myModel$results)) return(0)

	.womSetUseMH(myModel$settings$useMH)
	theWidget("analysisMHButton")$setActive(myModel$settings$useMH)
	if(myModel$settings$useMH)
	{
		theWidget("analysisMHScaleEntry")$setText(myModel$settings$MHScale)
		theWidget("analysisMHThinEntry")$setText(myModel$settings$MHThin)
	}else{
		theWidget("analysisLowerEpsEntry")$setText(myModel$settings$epsLow)
		theWidget("analysisUpperEpsEntry")$setText(myModel$settings$epsUpp)
		theWidget("analysisLeapfrogEntry")$setText(myModel$settings$leapfrog)
	}
	
	theWidget("analysisItersEntry")$setText(myModel$settings$MCMCIters)
	theWidget("analysisBurninEntry")$setText(myModel$settings$burninIters)
	theWidget("analysisStartOptimEntry")$setText(myModel$settings$startIters)
	theWidget("analysisPredProbButton")$setActive(myModel$settings$predProbs)
	theWidget("analysisAlarmButton")$setActive(myModel$settings$alarm)
}

.toggled_analysis_useMH <- function(toggle)
{
	checkbox = theWidget("analysisMHButton")
	toggled = gtkToggleButtonGetActive(checkbox)
	
	if(toggled){
		.womSetUseMH(TRUE)
	}else{
		.womSetUseMH(FALSE)
	}
}

.womSetUseMH <-function(useMH=FALSE)
{
	theWidget("analysisLowerEpsEntry")$setSensitive(!useMH)
	theWidget("analysisUpperEpsEntry")$setSensitive(!useMH)
	theWidget("analysisLeapfrogEntry")$setSensitive(!useMH)
	
	theWidget("analysisMHScaleEntry")$setSensitive(useMH)
	theWidget("analysisMHThinEntry")$setSensitive(useMH)
}

.clicked_analysis_analyze_selected<-function(button)
{
	allGood = .womCheckAnalysisSettings()
	if(!allGood) return(0)
	
	anythingToDo = .womCheckToBeAnalyzed()
	if(!anythingToDo) return(0)

	settings = .womCompileMCMCSettings()
	.womStartAnalyses(settings)
	
}

.womStartAnalyses<-function(settings)
{
	freezeGUI()
	on.exit(thawGUI())
	treeview <- theWidget("analysisDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	
	toBeAnalyzed = which(.womTreeModelGetNthCol(model,n=.womDefinedModelsTreeCols("toBeAnalyzed"))==1)
	numAnalyses = length(toBeAnalyzed)
	
	allProgressBar <- theWidget("analysisAllProgressbar")
	
	for(i in toBeAnalyzed)
	{
		index = which(i==toBeAnalyzed)
		gtkProgressBarSetFraction(allProgressBar,(index-1)/numAnalyses)
		wommbatAnalysis$Models[[i]]$settings <- settings
		.womAnalyzeOneModel(i)
		.womSetStatusBarText("Done.")
		wommbatAnalysis$Models[[i]]$results$timeAnalyzed = Sys.time()
		.womSetResultsInTableview(i)
	}
	gtkProgressBarSetFraction(allProgressBar,1)
	

	if(settings$alarm) alarm()
}

.womSetResultsInTableview	<-	function(modelNumber)
{
	treeview <- theWidget("modelsDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	
	if(is.character(modelNumber))
		modelNumber = which(names(wommbatAnalysis$Models)==modelNumber)
	
	.womActiveDiagnosticsTab(status=TRUE)
	.womActiveResultsTab(status=TRUE)
	.womCreateDiagnosticsDefinedModelsColumns()
	.womCreateResultsDefinedModelsColumns()
	
	myModel=wommbatAnalysis$Models[[modelNumber]]
	
	foreCol = .womChooseAccCol(myModel$results$accRate)
	
		
	iter = gtkTreeModelIterNthChild(model,parent=NULL, modelNumber-1)$iter
	model$set(iter, .womDefinedModelsTreeCols("hasResults"), TRUE,
					.womDefinedModelsTreeCols("toBeAnalyzed"), FALSE,
					.womDefinedModelsTreeCols("iterations"), as.integer(myModel$settings$MCMCIters),
					.womDefinedModelsTreeCols("effectiveIterations"), as.integer(myModel$settings$effectiveIters),
					.womDefinedModelsTreeCols("timeAnalyzed"), as.character(myModel$results$timeAnalyzed),
					.womDefinedModelsTreeCols("accRate"), as.numeric(.womSaneNum(myModel$results$accRate,2)),
					.womDefinedModelsTreeCols("useMet"), myModel$settings$useMH,
					.womDefinedModelsTreeCols("burnin"), as.integer(myModel$settings$burninIters),
					.womDefinedModelsTreeCols("DIC"), 
					as.numeric(.womSaneNum(myModel$results$DIC,1)),
					.womDefinedModelsTreeCols("DICtext"), 
					.womSaneNum(myModel$results$DIC,1),
					.womDefinedModelsTreeCols("pD"), as.numeric(.womSaneNum(myModel$results$pD,1)),
					.womDefinedModelsTreeCols("accColor"), foreCol,
					.womDefinedModelsTreeCols("logLikePostMean"),  as.numeric(.womSaneNum(myModel$results$logLikePostMean,1))
				)
	if(myModel$settings$useMH)
	{
		model$set(iter,
					.womDefinedModelsTreeCols("metropScale"), myModel$settings$MHScale,
					.womDefinedModelsTreeCols("metropThin"), myModel$settings$metropThin
				)
	}else
	{
		model$set(iter,
					.womDefinedModelsTreeCols("epsUpp"), as.character(myModel$settings$epsUpp),
					.womDefinedModelsTreeCols("epsLow"), as.character(myModel$settings$epsLow),
					.womDefinedModelsTreeCols("leapfrog"), as.integer(myModel$settings$leapfrog)
				)
	}
}

.womCheckAnalysisSettings<-function()
{
	startIters = try(as.integer(theWidget("analysisStartOptimEntry")$getText()),silent=TRUE)
	MCMCIters = try(as.integer(theWidget("analysisItersEntry")$getText()),silent=TRUE)
	burninIters = try(as.integer(theWidget("analysisBurninEntry")$getText()),silent=TRUE)

	if(inherits(startIters, "try-error") | is.na(startIters))
	{
		gWidgets::gmessage(paste("Could not coerce the optim() iterations to integer."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if(inherits(MCMCIters, "try-error") | inherits(burninIters, "try-error") |
		is.na(MCMCIters) | is.na(burninIters))
	{
		gWidgets::gmessage(paste("Could not coerce MCMC or burnin iterations to integers."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if(MCMCIters<1 | burninIters<1 | MCMCIters<=burninIters)
	{
		gWidgets::gmessage(paste("MCMC and burnin iterations must both be greater than 0, and burnin iterations must be less than the total MCMC iterations."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}

	MHcheckbox = theWidget("analysisMHButton")
	useMH = gtkToggleButtonGetActive(MHcheckbox)
	
	if(useMH)
	{
		MHScale = try(as.numeric(theWidget("analysisMHScaleEntry")$getText()),silent=TRUE)
		MHThin = try(as.integer(theWidget("analysisMHThinEntry")$getText()),silent=TRUE)
		
		if(inherits(MHScale, "try-error") | inherits(MHThin, "try-error") | 
			is.na(MHScale) | is.na(MHThin))
		{
			gWidgets::gmessage(paste("Could not coerce random walk Metropolis-Hastings settings to numerics."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(FALSE)
		}
		if(MHScale <= 0 | MHThin<1 | MHThin>=MCMCIters)
		{
			gWidgets::gmessage(paste("Metropolis-Hastings scale must be greater than 0, thinning must be at least 1 (but less than the number of total MCMC iterations)."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(FALSE)
		}
	
	}else
	{
		epsLow = try(as.numeric(theWidget("analysisLowerEpsEntry")$getText()),silent=TRUE)
		epsUpp = try(as.numeric(theWidget("analysisUpperEpsEntry")$getText()),silent=TRUE)
		leapfrog = try(as.integer(theWidget("analysisLeapfrogEntry")$getText()),silent=TRUE)
		if(inherits(epsLow, "try-error") | inherits(epsUpp, "try-error") | inherits(leapfrog, "try-error") |
			is.na(epsLow) | is.na(epsUpp) | is.na(leapfrog))
		{
			gWidgets::gmessage(paste("Could not coerce hybrid MCMC settings to numerics."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(FALSE)
		}
		if(epsLow<=0 | epsUpp<=0 | epsLow>epsUpp | leapfrog<1)
		{
			gWidgets::gmessage(paste("Lower and upper epsilon must be greater than 0, lower epsilon must be less than upper epsilon, and leapfrog steps must be at least 1."), title="Settings error",
			icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(FALSE)
		}
	
	}
	return(TRUE)
}

.womCheckToBeAnalyzed<-function()
{
	treeview <- theWidget("analysisDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	
	toBeAnalyzed = .womTreeModelGetNthCol(model,n=.womDefinedModelsTreeCols("toBeAnalyzed"))
	if(any(toBeAnalyzed==1)){
		.womSetStatusBarText("Starting analysis...")
		return(TRUE)
	}else{
		.womSetStatusBarText("Nothing to do.")
		return(FALSE)
	}
	
}

.womCompileMCMCSettings <- function()
{
	settings = list()
	
	# Entries
	settings$startIters = as.integer(theWidget("analysisStartOptimEntry")$getText())
	settings$MCMCIters = as.integer(theWidget("analysisItersEntry")$getText())
	settings$burninIters = as.integer(theWidget("analysisBurninEntry")$getText())
	
	settings$epsLow = as.numeric(theWidget("analysisLowerEpsEntry")$getText())
	settings$epsUpp = as.numeric(theWidget("analysisUpperEpsEntry")$getText())
	settings$leapfrog = as.integer(theWidget("analysisLeapfrogEntry")$getText())

	settings$MHScale = as.numeric(theWidget("analysisMHScaleEntry")$getText())
	settings$MHThin = as.integer(theWidget("analysisMHThinEntry")$getText())

	# toggles
	MHcheckbox = theWidget("analysisMHButton")
	settings$useMH = gtkToggleButtonGetActive(MHcheckbox)
	predCheckbox = theWidget("analysisPredProbButton")
	settings$predProbs = gtkToggleButtonGetActive(predCheckbox)
	alarmCheckbox = theWidget("analysisAlarmButton")
	settings$alarm = gtkToggleButtonGetActive(alarmCheckbox)

	return(settings)
}

.womFormatOutput<-function(output,modelNum)
{
	myModel = wommbatAnalysis$Models[[modelNum]]
	cov = myModel$model$covNgroups>0
	results=list()

	if(cov)
	{
		results$meanChains = coda::mcmc(t(output[[3]]))
		results$covChains = list()
        results$corChains = list()
        if(myModel$settings$predProbs) results$predVals = output[[5]]
        .womSetStatusBarText("Formatting cov/cor matrices...")
        for(i in 1:myModel$model$covNgroups)
          {
            results$covChains[[i]]  = apply(output[[4]][[i]],3,solve)
			dim(results$covChains[[i]]) = c(myModel$model$covSizes[i],myModel$model$covSizes[i],as.integer(myModel$settings$effectiveIters))
			results$corChains[[i]]  = apply(results$covChains[[i]],3,cov2cor)
			dim(results$corChains[[i]]) = dim(results$covChains[[i]])
		  }
		
		
	}else
	{
		if(myModel$settings$predProbs) results$predVals = output[[3]] 
	}
	
	results$accRate = round(mean(diff(output[[2]])!=0),3)
	postSD = apply(output[[1]][,(myModel$settings$burninIters+1):as.integer(myModel$settings$effectiveIters)],1,sd)
	effSize = apply(output[[1]][,(myModel$settings$burninIters+1):as.integer(myModel$settings$effectiveIters)],1,coda::effectiveSize)
	specar = unlist(apply(output[[1]][,(myModel$settings$burninIters+1):as.integer(myModel$settings$effectiveIters)],1,coda::spectrum0.ar))
	specar = matrix(specar,nrow=2)[1,]
	results$pointEst = .womNiceParVec(rowMeans(output[[1]][,(myModel$settings$burninIters+1):as.integer(myModel$settings$effectiveIters)]), 
										myModel$model$newDat2Cat, myModel$model$newDat2Cont, myModel$model$namedDat2, myModel$model$effects, 
										myModel$model$incCont, TRUE,	
										postSD,sqrt(specar/(myModel$settings$effectiveIters-myModel$settings$burninIters)))
	chainnames=paste(results$pointEst[,4],results$pointEst[,2],results$pointEst[,3],sep=" ")
	chainnames=paste(chainnames,results$pointEst[,1],sep=" on ")
	## Add in names of Covariance matrices!
	results$effectChains=coda::mcmc(t(output[[1]]))
	dimnames(results$effectChains)[[2]] = chainnames
	if(myModel$settings$useMH){
		results$likeChain = output[[2]][((1:as.integer(myModel$settings$MCMCIters))%%as.integer(myModel$settings$MHThin))==0]
	}else{
		results$likeChain = output[[2]]
	}
	likePostMean <- womRlogLikelihood(results$pointEst[,5],myModel)
	results$logLikePostMean <- likePostMean
	.womSetStatusBarText("Computing DIC...")
	DICpD = .womComputeDIC(results$likeChain[(myModel$settings$burninIters+1):as.integer(myModel$settings$effectiveIters)],likePostMean)
	results$DIC = DICpD[1]
	results$pD = DICpD[2]
	
	if(cov){
		covNames <- paste(results$pointEst[,2],results$pointEst[,1],sep=" on ")[myModel$model$parStart+1]
		dim(covNames) <- dim(myModel$model$parStart)
		for(i in 1:myModel$model$covNgroups)
		{
			dimnames(results$covChains[[i]])[[1]] = covNames[i,1:myModel$model$covSizes[i]]
			dimnames(results$corChains[[i]])[[1]] = covNames[i,1:myModel$model$covSizes[i]]
			dimnames(results$covChains[[i]])[[2]] = covNames[i,1:myModel$model$covSizes[i]]
			dimnames(results$corChains[[i]])[[2]] = covNames[i,1:myModel$model$covSizes[i]]

		}
	}	
	return(results)
}


.womAnalyzeOneModel <- function(modelNumber)
{
	myModel <- wommbatAnalysis$Models[[modelNumber]]
	wommbatAnalysis$Models[[modelNumber]]$messages = list()	

	
	startingvals = jitter((1:(sum(myModel$model$effects)+3))*0)
    startingvals[1] = myModel$priors$muKMean
    startingvals[2] = myModel$priors$muZMean
    startingvals[3] = myModel$priors$muGMean
		
	optimIters <- myModel$settings$startIters
	if(optimIters>0)
	{
		.womSetStatusBarText("Using optim() to find starting values...")
		optimOut = optim(startingvals,womRlogPosterior,womRgradLogPosterior,method="BFGS",control=list(maxit=optimIters),setup=myModel,hessian=TRUE) 
	   	wommbatAnalysis$Models[[modelNumber]]$messages$optimConvergence = optimOut$convergence
		startingvals=optimOut$par
		.womSetStatusBarText("Inverting Hessian...")
		wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights = diag(solve(optimOut$hessian))
		if(any(wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights<=0)){
			wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights = wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights*0 + 1
		}
	}else
	{
		wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights = startingvals*0 + 1	
	}	
	
	if(myModel$settings$useMH){
		wommbatAnalysis$Models[[modelNumber]]$settings$metropSD = 
		sqrt(wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights)*myModel$settings$MHScale
		wommbatAnalysis$Models[[modelNumber]]$settings$effectiveIters = floor(myModel$settings$MCMCIters/myModel$settings$MHThin)
	}else{
		wommbatAnalysis$Models[[modelNumber]]$settings$metropSD = 0				
		wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights = wommbatAnalysis$Models[[modelNumber]]$settings$MCMCweights*0 + 1
		#pack@settings@MCMCweights = pack@settings@MCMCweights/exp(mean(log(pack@settings@MCMCweights)))
		wommbatAnalysis$Models[[modelNumber]]$settings$effectiveIters = myModel$settings$MCMCIters
	}
    
	wommbatAnalysis$Models[[modelNumber]]$settings$startingValues <- startingvals
	
	.womSetStatusBarText("Generating MCMC chains.")
	if(myModel$model$covNgroups>0){
		.womAnalysisCallWithCov(modelNumber)
	}else{
		.womAnalysisCallNoCov(modelNumber)
	}
	
}

.womAnalysisCallWithCov <- function(modelNum)
{
	myModel=wommbatAnalysis$Models[[modelNum]]
	
	singleProgressBar <- theWidget("analysisCurrentProgressbar")
	gtkProgressBarSetFraction(singleProgressBar,0)
	
	progress = ceiling(myModel$settings$MCMCIters/100)
	pbFun = function(samps)
	{
		gtkProgressBarSetFraction(singleProgressBar,samps/myModel$settings$MCMCIters)
	}
	output = .Call("WM2_GibbsSampler", 
				as.integer(myModel$settings$MCMCIters), 
				myModel$settings$startingValues, 
				as.integer(myModel$model$covNgroups),
				as.integer(myModel$model$covObs), 
				as.integer(myModel$model$covSizes),
				as.integer(myModel$model$parStart),
				as.integer(myModel$model$effSlope),
				as.integer(myModel$model$newDat2Cat[,1]), 
				as.integer(myModel$model$newDat2Cat[,2]), 
				as.integer(myModel$model$newDat2Cat[,3]), 
				as.integer(myModel$model$newDat2Cat[,4]), 
				as.integer(as.character(myModel$model$newDat2Cat[,5])), 
				.womIntegerMatrix(myModel$model$newDat2Cat[,-(1:5)]), 
				as.matrix(myModel$model$newDat2Cont[,-(1:5)]), 
				as.integer(myModel$model$effects[1,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$inCovMat[1,]), 
				as.integer(myModel$model$effects[2,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$inCovMat[2,]), 
				as.integer(myModel$model$effects[3,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$inCovMat[3,]), 
				myModel$priors$IGa0, myModel$priors$IGb0, 
				myModel$priors$muKMean, myModel$priors$muKSD^2, 
				myModel$priors$muZMean, myModel$priors$muZSD^2, 
				myModel$priors$muGMean, myModel$priors$muGSD^2, 
				1, as.integer(myModel$Ktype), myModel$settings$epsLow, 
				myModel$settings$epsUpp - myModel$settings$epsLow, 
				as.integer(myModel$settings$leapfrog), 
				myModel$settings$MCMCweights, myModel$priors$invWishartScalar,
				as.integer(progress),pbFun,new.env(),
				as.integer(myModel$settings$predProbs),as.integer(myModel$settings$useMH),
				myModel$settings$metropSD,as.integer(myModel$settings$MHThin),
				package = "WMCapacity")
				
	gtkProgressBarSetFraction(singleProgressBar,1)	  
	#wommbatAnalysis$Models[[modelNum]]$rawOutput = output
	wommbatAnalysis$Models[[modelNum]]$results=.womFormatOutput(output,modelNum)
	
}		  

.womAnalysisCallNoCov <- function(modelNum)
{
	myModel=wommbatAnalysis$Models[[modelNum]]
	
	singleProgressBar <- theWidget("analysisCurrentProgressbar")
	gtkProgressBarSetFraction(singleProgressBar,0)
	
	progress = ceiling(myModel$settings$MCMCIters/100)
	pbFun = function(samps)
	{
		gtkProgressBarSetFraction(singleProgressBar,samps/myModel$settings$MCMCIters)
	}
	output = .Call("WM2_GibbsSamplerNoCov", 
				as.integer(myModel$settings$MCMCIters), 
				myModel$settings$startingValues, 
				as.integer(myModel$model$newDat2Cat[,1]), 
				as.integer(myModel$model$newDat2Cat[,2]), 
				as.integer(myModel$model$newDat2Cat[,3]), 
				as.integer(myModel$model$newDat2Cat[,4]), 
				as.integer(as.character(myModel$model$newDat2Cat[,5])), 
				.womIntegerMatrix(myModel$model$newDat2Cat[,-(1:5)]), 
				as.matrix(myModel$model$newDat2Cont[,-(1:5)]), 
				as.integer(myModel$model$effects[1,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$effects[1,]*0), 
				as.integer(myModel$model$effects[2,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$effects[2,]*0), 
				as.integer(myModel$model$effects[3,]), 
				as.integer(myModel$model$incCont), 
				as.integer(myModel$model$effects[3,]*0), 
				myModel$priors$IGa0, myModel$priors$IGb0, 
				myModel$priors$muKMean, myModel$priors$muKSD^2, 
				myModel$priors$muZMean, myModel$priors$muZSD^2, 
				myModel$priors$muGMean, myModel$priors$muGSD^2, 
				1, as.integer(myModel$Ktype), myModel$settings$epsLow, 
				myModel$settings$epsUpp - myModel$settings$epsLow, 
				as.integer(myModel$settings$leapfrog), 
				myModel$settings$MCMCweights,
				as.integer(progress),pbFun,new.env(),
				as.integer(myModel$settings$predProbs),as.integer(myModel$settings$useMH),
				myModel$settings$metropSD,as.integer(myModel$settings$MHThin),
				package = "WMCapacity")
	gtkProgressBarSetFraction(singleProgressBar,1)	  
	#wommbatAnalysis$Models[[modelNum]]$rawOutput = output
	wommbatAnalysis$Models[[modelNum]]$results=.womFormatOutput(output,modelNum)
}


.womNiceParVec=function(par,newDat2Cat,newDat2Cont,namedDat2,effects,incCont,useA,sd=NA,MCMCerr=NA)
{  
  newPars=data.frame(array(par,dim=c(length(par),7)))
  colnames(newPars)=c("Parameter","Effect","Level","Type","Estimate","Post. SD","Estimate MCMC Error")
  newPars[,5]=par
  newPars[,6]=sd
  newPars[,7]=MCMCerr
  newPars[1,1:4]=c("k","grand mean","-","-")
  intcolnames1=colnames(newDat2Cat)[-(1:5)]
  intcolnames2=colnames(newDat2Cont)[-(1:5)]
  intcolnames=paste(intcolnames1,intcolnames2,sep=".x.")
  intcolnames[incCont==0]=intcolnames1[incCont==0]
  parnames=c("k","z","g")
  
  if(useA){
    newPars[2,1:4]=c("z","grand mean","-","-")
    newPars[3,1:4]=c("g","grand mean","-","-")
    counter=4
    for(i in 1:dim(effects)[2])
  	for(j in 1:3){
	      if(effects[j,i]>0){
		for(k in 1:effects[j,i]){
		      cellname=namedDat2[(k==newDat2Cat[,i+5]),i+5][1]
		      newPars[counter,1:4]=c(parnames[j],intcolnames[i],as.character(cellname),c("Random Effect","Slope")[incCont[i]+1])
		      counter=counter+1
		}
	      }
		      	      
	}
  }else{
    newPars[2,1:4]=c("g","grand mean","-","-")
    counter=3
       for(i in 1:dim(effects)[2])
  	for(j in 1:3){
	      if(effects[j,i]>0){
		for(k in 1:effects[j,i]){
		      cellname=namedDat2[(k==newDat2Cat[,i+5]),i+5][1]
		      newPars[counter,1:4]=c(parnames[j],intcolnames[i],as.character(cellname),c("Random Effect","Slope")[incCont[i]+1])
		      counter=counter+1
		}
	      }
		      	      
	}
  }
newPars
}

