#
# Author: Ian Fellows
###############################################################################

.factorAnalysisCheckFunction <- function(state){
	#make sure at least two variables are selected
	if(length(state$variables)<2)
		return("Please select at least two variables")
	return("")
}

.factorAnalysisRunFunction <- function(state){
	#print(state) #a print statement is useful for debugging
	
	#make formula
	form <-paste( " ~ " , state$variables[1])
	for( var in state$variables[-1])
		form <- paste(form,"+",var)
	
	#make prcomp call
	cmd <- paste("pr.model <-prcomp(", form, ",", state$data)
	if("Center" %in%state$Transformation)
		cmd <- paste(cmd,", center=TRUE")
	if("Scale" %in%state$Transformation)
		cmd <- paste(cmd,",scale=TRUE")
	cmd <- paste(cmd,")")
	
	#always print model
	cmd <- paste (cmd,"\n","print(pr.model)")
	
	#output summary and plot if asked for
	if("Summary" %in% state$Output)
		cmd <- paste(cmd,"\n","summary(pr.model)")
	if("Scree Plot" %in% state$Output)
		cmd <- paste(cmd,"\n","screeplot(pr.model)")
	
	#execute command as if typed into console
	execute(cmd)
}

makeFactorAnalysisDialog <- function(){
	#make dialog
	dialog <- new(SimpleRDialog)
	dialog$setSize(500L,400L)
	dialog$setTitle("Factor Analysis")
	
	#add variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(dialog,variableSelector,10,400,850,10)
	
	#add a list for the variables
	variableList<- new(VariableListWidget,variableSelector)
	variableList$setTitle("variables")
	addComponent(dialog, variableList,100,900,450, 420)
	
	#options for transforming the variables
	transBoxes <- new(CheckBoxesWidget,"Transformation",c("Center","Scale"))
	addComponent(dialog, transBoxes,500,900,670, 540)
	transBoxes$setDefaultModel(c("Scale"))
	
	#output options
	outBoxes <- new(CheckBoxesWidget,"Output",c("Summary","Scree Plot"))
	addComponent(dialog, outBoxes,680,900,850, 540)
	dialog$setCheckFunction(toJava(.factorAnalysisCheckFunction))
	dialog$setRunFunction(toJava(.factorAnalysisRunFunction))
	return(dialog)
}



makeListExampleDialog <-function(){
	leftList <- new(ListWidget,"Left hand list",c("Item 1","Item 2","Item 3","Item 4"))
	setSize(leftList,150,210)
	rightList <- new(ListWidget,"Right hand list")
	setSize(rightList,150,210)
	buttons <- new(AddRemoveButtons,leftList,rightList)
	
	dialog <- new(SimpleRDialog)
	dialog$setSize(500L,300L)
	addComponent(dialog,leftList,10,400,600,450,,,"NONE","NONE")
	addComponent(dialog,rightList,10,900,600,450,,,"NONE","NONE")
	addComponent(dialog,buttons,100,600,600,450,"NONE","NONE")
	
	runDialog <- function(state) print(state)
	dialog$setRunFunction(toJava(runDialog))
	
	return(dialog)
}


makeScatterPlotDialog <- function(){
	PlotRDialog <- J("example.PlotRDialog")
	.scatterPlotDialog <- new(PlotRDialog)
	return(.scatterPlotDialog)
}

.onLoad <- function(libname, pkgname) { 
	
	#if deducer gui is not running, do minimal load
	deducerLoaded <- try(.deducer == .jnull(),silent=TRUE)
	if(inherits(deducerLoaded,"try-error") || deducerLoaded)
		return(NULL)
	
	#loads example.jar
	.jpackage(pkgname,lib.loc=libname)
	
	.registerDialog("Factor Analysis", makeFactorAnalysisDialog)
	.registerDialog("Scatter Plot", makeScatterPlotDialog)
	.registerDialog("List Example", makeListExampleDialog)
	
	#add menu items
	deducer.addMenu("Example")
	deducer.addMenuItem("Factor Analysis",,".getDialog('Factor Analysis')$run()","Example")
	deducer.addMenuItem("Scatter Plot",,".getDialog('Scatter Plot')$run()","Example")
	deducer.addMenuItem("Widget display",,"J('example.ExampleDialog')$run()","Example")
	deducer.addMenuItem("List Example",,".getDialog('List Example')$run()","Example")
	if(.windowsGUI){
		winMenuAdd("Example")
		winMenuAddItem("Example", "Factor Analysis", "deducer('Factor Analysis')")
		winMenuAddItem("Example", "Scatter Plot", "deducer('Scatter Plot')")
		winMenuAddItem("Example", "Widget Display", "deducer('Widget display')")
		winMenuAddItem("Example", "List Example", "deducer('List Example')")
	}else if(.jgr){
		jgr.addMenu("Example")
		jgr.addMenuItem("Example", "Factor Analysis", "deducer('Factor Analysis')")
		jgr.addMenuItem("Example", "Scatter Plot", "deducer('Scatter Plot')")
		jgr.addMenuItem("Example", "Widget Display", "deducer('Widget display')")
		jgr.addMenuItem("Example", "List Example", "deducer('List Example')")
	}
}