.womSetupDiagnosticsComboBox<-function()
{
	itersSpace = theWidget('diagnosticItersComboSpace')
	typeSpace = theWidget('diagnosticTypeComboSpace')
	StateEnv$itersCombo = gtkComboBoxNewText()
	StateEnv$typeCombo = gtkComboBoxNewText()
	gtkComboBoxAppendText(StateEnv$itersCombo, "All")
	gtkComboBoxAppendText(StateEnv$typeCombo, "effects")
	itersSpace$packStart(StateEnv$itersCombo,FALSE,FALSE,0)
	typeSpace$packStart(StateEnv$typeCombo,FALSE,FALSE,0)
	gtkComboBoxSetActive(StateEnv$itersCombo,0)
	gtkComboBoxSetActive(StateEnv$typeCombo,0)
	#Connect chain limit signal
	StateEnv$handlers$diagnosticLimitCombo <- gSignalConnect(StateEnv$itersCombo, "changed", .selected_number_chain_iterations)
	# Connect parameter type signal
	StateEnv$handlers$diagnosticTypeCombo <- gSignalConnect(StateEnv$typeCombo, "changed", .selected_diagnostic_parameter_type)
}


.womActiveDiagnosticsTab<-function(status=TRUE)
{

	theWidget("diagnosticPageBox")$setSensitive(status)

}


.womCreateDiagnosticsDefinedModelsColumns <- function()
{
	#modelsTreeview <- theWidget("modelsDefinedModelsTreeview")
	diagnosticsTreeview <- theWidget("diagnosticDefinedModelsTreeview")
	
	#model = gtkTreeViewGetModel(modelsTreeview)	
	#gtkTreeViewSetModel(diagnosticsTreeview, model)	
	
	
	treeview = diagnosticsTreeview
	
	#Create tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		
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
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("accRate"))
		
		
		# iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("iterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("iterations"))
		
		# effective iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effective Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("effectiveIterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("effectiveIterations"))
		
		# burnin
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Burnin", renderer, 
  								text = .womDefinedModelsTreeCols("burnin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("burnin"))
		
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

.womGetDiagnosticsModelSelection<-function()
{
	
	treeview = theWidget("diagnosticDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	
		
	# Next line caused crash - "visible" is not a valid column!
	# I don't remember why I used this. (cut and paste?)
	# vis <- model$get(iter, .womDefinedModelsTreeCols("visible"))[[1]]
	vis <- NULL
	
	return(list(iter=iter,visible=vis))
}


.selected_diagnostics_model_row <- function(treeview, path, view_column)
{
		
	treeview <- theWidget("diagnosticDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetDiagnosticsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	
	myModel = wommbatAnalysis$Models[[modelName]]

	limitCombo = StateEnv$itersCombo
	typeCombo = StateEnv$typeCombo

	
	if(is.null(myModel$results))
	{
		.womClearDiagnosticPlots()
		limitCombo$setSensitive(FALSE)
		typeCombo$setSensitive(FALSE)
		theWidget("hscrollbar1")$setSensitive(FALSE)
		return()
	}
	limitCombo$setSensitive(TRUE)
	typeCombo$setSensitive(TRUE)
	theWidget("hscrollbar1")$setSensitive(TRUE)
	
	
	# Combo box for iterations limit
	currentValue = limitCombo$getActiveText()
	currentValue = ifelse(is.null(currentValue),"All",currentValue)
	
	gSignalHandlerBlock(limitCombo, StateEnv$handlers$diagnosticLimitCombo)
	burnin <- myModel$settings$burninIters
	upperLimit <- myModel$settings$effectiveIters
	

	# Make options model
	limitOptions <-c(upperLimit-burnin,500,1000,5000)
	limitOptions <- limitOptions[limitOptions <= upperLimit-burnin]
	limitOptions[1]="All"
	findCurrent = match(currentValue,limitOptions)
	selectIters = ifelse(is.na(findCurrent),0,findCurrent-1)
	
	clearComboModel(StateEnv$itersCombo)
	
	for (LO in limitOptions) {
		gtkComboBoxAppendText(limitCombo,LO)
	}
	gtkComboBoxSetActive(limitCombo,selectIters)
		
	# Combo box for parameter type
	#currentValueType = theWidget("diagnosticTypeComboBox")$getActiveText()
	currentValueType = typeCombo$getActiveText()
	currentValueType = ifelse(is.null(currentValueType),"effect",currentValueType)
	selectType=0
	
	gSignalHandlerBlock(typeCombo, StateEnv$handlers$diagnosticTypeCombo)
	
	clearComboModel(StateEnv$typeCombo)
	gtkComboBoxAppendText(typeCombo,"effect")
	
	if(myModel$model$covNgroups>0){
		gtkComboBoxAppendText(typeCombo,"correlation")
		if(currentValueType=="correlation") selectType=1
	}
	gtkComboBoxSetActive(typeCombo,selectType)

	
	maxPar <- length(myModel$results$pointEst[,1])
	
	scrollBar <- theWidget("hscrollbar1")
	gSignalHandlerBlock(scrollBar, StateEnv$handlers$diagnosticScrollBar1)
	adjustment <- gtkAdjustmentNew(value = 1, lower = 1, upper = maxPar, step.incr=1)
	scrollBar$setAdjustment(adjustment)
	gtkAdjustmentSetValue(scrollBar$getAdjustment(), 1)
			
	.womDisplayNthDiagnosticPlot(1)
	
	gSignalHandlerUnblock(scrollBar, StateEnv$handlers$diagnosticScrollBar1)
	gSignalHandlerUnblock(typeCombo,StateEnv$handlers$diagnosticTypeCombo)
	gSignalHandlerUnblock(limitCombo, StateEnv$handlers$diagnosticLimitCombo)
}

.womClearDiagnosticPlots <- function()
{
	if(!is.null(StateEnv$Graphics$ChainsDevice)){
		dev.set(which = StateEnv$Graphics$ChainsDevice)
		plot.new()
	}
	if(!is.null(StateEnv$Graphics$DensityDevice)){
		dev.set(which = StateEnv$Graphics$DensityDevice)
		plot.new()
	}
	if(!is.null(StateEnv$Graphics$ACFDevice)){
		dev.set(which = StateEnv$Graphics$ACFDevice)
		plot.new()
	}
}

.womGetApproprateChain <- function()
{
	treeview <- theWidget("diagnosticDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetDiagnosticsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	myModel = wommbatAnalysis$Models[[modelName]]
	
	scrollBar <- theWidget("hscrollbar1")
	n = floor(scrollBar$getValue())
	if(n==0) n=1
	
	parameterType <- StateEnv$typeCombo#theWidget("diagnosticTypeComboBox")
	typeModel <- parameterType$getModel()
	
	typeIdx <- gtkComboBoxGetActiveIter(parameterType)
	if(!typeIdx$retval)
	{
		myType <- "effects"
	}else
	{
		iter <- typeIdx$iter
		myType <- gtkTreeModelGetValue(typeModel,iter,0)$value
	}
	
	if(myType=="correlation")
	{
		if(myModel$model$covNgroups==0) return(NULL)
	
		nParameters <- myModel$model$covSizes*(myModel$model$covSizes-1)/2
	
		if(n>sum(nParameters)) return(NULL)
	
		group <- which(n<=cumsum(nParameters))[1]
		idxInGroup = ifelse(group==1,n,n-cumsum(nParameters)[group-1])
		groupPars = nParameters[group]
		groupSize <- myModel$model$covSizes[group]
	
		x=matrix(1:groupSize^2,groupSize,groupSize)
		idx = x[lower.tri(x)][idxInGroup]
		row = (idx-1)%%groupSize + 1
		column = (idx-1)%/%groupSize + 1
		lab = paste("Correlation: ",dimnames(myModel$results$corChains[[group]])[[1]][row]," and ",
					dimnames(myModel$results$corChains[[group]])[[1]][column],sep="")
		
		return(list(chain=myModel$results$corChains[[group]][row,column,],group=group,row=row,column=column,label=lab))
	}else 
	{
		if(n>length(myModel$results$effectChains[1,])) return(0)
		
		parameterName <- paste(myModel$results$pointEst[n,2]," on ",myModel$results$pointEst[n,1],": ",myModel$results$pointEst[n,3],sep="")
		chain <- myModel$results$effectChains[,n]
		return(list(chain=chain,group=NULL,row=NULL,column=NULL,label=parameterName))
	}
	return(NULL)	
}

.womDisplayNthDiagnosticPlot <- function(n)
{
	treeview <- theWidget("diagnosticDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetDiagnosticsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	
	myModel = wommbatAnalysis$Models[[modelName]]
	parameterLabel <- theWidget("diagnosticsParameterLabel")

	parameterType <- StateEnv$typeCombo#theWidget("diagnosticTypeComboBox")
	typeModel <- parameterType$getModel()
	
	typeIdx <- gtkComboBoxGetActiveIter(parameterType)
	if(!typeIdx$retval)
	{
		myType <- "effects"
	}else
	{
		iter <- typeIdx$iter
		myType <- gtkTreeModelGetValue(typeModel,iter,0)$value
	}
	
	upperLimit <- myModel$settings$MCMCIters
	
	burnin <- myModel$settings$burninIters
	chainData <- .womGetApproprateChain()

	parameterName <- chainData$lab
	gtkLabelSetText(parameterLabel, parameterName)
	
	
	chain <- chainData$chain
	postMean <- mean(chain[-(1:burnin)])
	postQuant <- quantile(chain[-(1:burnin)],p=c(.025,.975))
	
	#Chain
	limitCombo <- StateEnv$itersCombo#theWidget("diagnosticCombobox1")
	limitModel <- gtkComboBoxGetModel(limitCombo)
	upperLimitIdx <- gtkComboBoxGetActiveIter(limitCombo)
	if(!upperLimitIdx$retval)
	{
		upperLimit <- myModel$settings$effectiveIters
	}else
	{
		iter <- upperLimitIdx$iter
		upperLimit <- gtkTreeModelGetValue(limitModel,iter,0)$value
		if(upperLimit=="All")
		{
			upperLimit <- myModel$settings$effectiveIters
		}else
		{
			upperLimit <- as.numeric(upperLimit) + burnin
		}
	}
	
	
	.womUpdateChainPlot(myModel,upperLimit)
	
	#Density
	if(!is.null(StateEnv$Graphics$DensityDevice)){
		dev.set(which = StateEnv$Graphics$DensityDevice)
	}else{
		da = theWidget("diagnosticDrawingarea2")
		cairoDevice::asCairoDevice(da)
		StateEnv$Graphics$DensityDevice = dev.cur()
	}
	par(mar=c(4.5,0,0,0))
	plot(density(chain[-(1:burnin)]),axes=FALSE,xlab="Parameter Estimate",main="")
	abline(v=postMean,col="blue",lwd=2)
	abline(v=postQuant,col="blue",lty=2)
	axis(1)
	box()
	
	#ACF
	if(!is.null(StateEnv$Graphics$ACFDevice)){
		dev.set(which = StateEnv$Graphics$ACFDevice)
	}else{
		da = theWidget("diagnosticDrawingarea3")
		cairoDevice::asCairoDevice(da)
		StateEnv$Graphics$ACFDevice = dev.cur()
	}
	par(mar=c(4.5,4,1,1))
	acf(chain[-(1:burnin)])
	axis(1)
	box()

		
}

.selected_diagnostic_parameter_type <- function(treeview, path)
{
	treeview <- theWidget("diagnosticDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetDiagnosticsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	myModel = wommbatAnalysis$Models[[modelName]]

	#parameterType <- theWidget("diagnosticTypeComboBox")
	parameterType <- StateEnv$typeCombo
	typeModel <- parameterType$getModel()
	
	typeIdx <- gtkComboBoxGetActiveIter(parameterType)
	if(!typeIdx$retval)
	{
		myType <- "effects"
	}else
	{
		iter <- typeIdx$iter
		myType <- gtkTreeModelGetValue(typeModel,iter,0)$value
	}
	
	if(myType=="correlation")
	{
		maxPar <- sum(myModel$model$covSizes*(myModel$model$covSizes-1)/2)
	}else
	{
		maxPar <- length(myModel$results$pointEst[,1])	
	}
	
	scrollBar <- theWidget("hscrollbar1")
	adjustment <- gtkAdjustmentNew(value = 1, lower = 1, upper = maxPar, step.incr=1)
	scrollBar$setAdjustment(adjustment)
	gtkAdjustmentSetValue(scrollBar$getAdjustment(), 1)
			
	.womDisplayNthDiagnosticPlot(1)	
}


.selected_number_chain_iterations <- function(treeview, path)
{
	
	treeview <- theWidget("diagnosticDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetDiagnosticsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	myModel = wommbatAnalysis$Models[[modelName]]
	
	burnin <- myModel$settings$burninIters
	
	#limitCombo <- theWidget("diagnosticCombobox1")
	limitCombo <- StateEnv$itersCombo
	limitModel <- gtkComboBoxGetModel(limitCombo)
	upperLimitIdx <- gtkComboBoxGetActiveIter(limitCombo)
	if(!upperLimitIdx$retval)
	{
		upperLimit <- myModel$settings$effectiveIters
	}else
	{
		iter <- upperLimitIdx$iter
		upperLimit <- gtkTreeModelGetValue(limitModel,iter,0)$value
		if(upperLimit=="All")
		{
			upperLimit <- myModel$settings$effectiveIters
		}else
		{
			upperLimit <- as.numeric(upperLimit) + burnin
		}
	}
	
	.womUpdateChainPlot(myModel,upperLimit)
}


.womUpdateChainPlot<-function(myModel,upperLimit)
{	
	scrollBar <- theWidget("hscrollbar1")
	n = floor(scrollBar$getValue())
	if(n==0) n=1
	
	burnin <- myModel$settings$burninIters
	
	chainData <- .womGetApproprateChain()
	chain <- chainData$chain
	postMean <- mean(chain[-(1:burnin)])
	postQuant <- quantile(chain[-(1:burnin)],p=c(.025,.975))
		
	# plots
	# Chains plot
	if(!is.null(StateEnv$Graphics$ChainsDevice)){
		dev.set(which = StateEnv$Graphics$ChainsDevice)
	}else{
		da = theWidget("diagnosticDrawingarea1")
		cairoDevice::asCairoDevice(da)
		StateEnv$Graphics$ChainsDevice = dev.cur()
	}
	par(mar=c(4.5,4,1,1))
	plot(chain[(burnin+1):upperLimit],axes=FALSE,ty='l',ylab="Parameter Estimate")
	abline(h=postMean,col="blue",lwd=2)
	abline(h=postQuant,col="blue",lty=2)
	lines(lowess(chain[(burnin+1):upperLimit],f=.25),col="red",lwd=2)
	axis(2)
	axis(1,at=c(1,upperLimit-burnin),labels=c(burnin+1,upperLimit))
	box()
	
}

.scrolled_diagnostics_scrollbar <- function(scrollBar,...)
{
	scrollBar <- theWidget("hscrollbar1")
	
	value = floor(scrollBar$getValue())
	.womDisplayNthDiagnosticPlot(value)
	
}

