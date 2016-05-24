
.womColumnsTreeCols <- function(name)
{
COLUMNS <- c(columnname=0,visible=1)
as.integer(COLUMNS[name])
}

.womColumnsOfInterestTreeCols <- function(name)
{
COLUMNS <- c(columnname=0,type=1,visible=2)
as.integer(COLUMNS[name])
}

.wom_clicked_open_csv <- function(button)
{
	#freezeGUI()
	#on.exit(thawGUI())
	fileChoose(action=".womOpenCSVFile", 
		type="open", 
		text="Select a CSV file...", 
		filter = list("CSV files" = list(patterns = c("*.csv","*.CSV")),
					  "text files" = list(mime.types = c("text/plain")),
					  "All files" = list(patterns = c("*")) 
					 )
		)
	
	#filename <- try(file.choose(), silent=T)
	#StateEnv$win$present()
	#if (inherits(filename, "try-error")) return()
	
	#.womOpenCSVFile(filename)
}

.womOpenCSVFile <- function(filename)
{
	dataset = try(read.csv(filename), silent=T)
	if (inherits(dataset, "try-error")){
		gWidgets::gmessage(paste("The file '",filename,"' could not be loaded. The probable cause this of error is that the file is not in CSV format.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	
	
	.womSetDataForColumnSelection(dataset,filename=filename)
	
}


.womCreateColumnSelectionColumns <- function()
{
	treeview1 <- theWidget("dataTreeview1")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	#Create left columns tree view
	if(is.null(gtkTreeViewGetColumn(treeview1,0))){
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview1$insertColumnWithAttributes(-1, "Column name", renderer, 
  								text = .womColumnsTreeCols("columnname") ,
								visible = .womColumnsTreeCols("visible"))
								
		column <- treeview1$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	}
}

.womCopyModelInfo <- function(project)
{
	rm(list=ls(envir=wommbatAnalysis),envir=wommbatAnalysis)
	x = ls(envir=project)
	eval(parse(text = paste("wommbatAnalysis$",x,"<-project$",x,sep="")))
}

.womLoadProject <- function(project, projectFile="")
{

	
	myPage <- gtkNotebookGetCurrentPage(theWidget("notebook1"))
	
	if(!is.environment(project)){
		gWidgets::gmessage(paste("The project could not be loaded. The object passed to the function is not an R environment.",
						sep=""), title="Project load error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)	
	}else{
		
		.womCopyModelInfo(project)
			
		# check data set
		if(is.null(wommbatAnalysis$data))
		{
			gWidgets::gmessage(paste("The project could not be loaded. No data was found; is the environment a WoMMBAT data analysis?",
				sep=""), title="Project load error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(0)	

		}
		
		# Load data
		.womSetDataForColumnSelection(wommbatAnalysis$data,removeModels=FALSE,filename=projectFile)
		gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("data"))
	
		if(!is.null(wommbatAnalysis$Ktype))
		{
			if(wommbatAnalysis$Ktype==0){
				theWidget("dataDesignPashlerRadio")$setActive(0)
				theWidget("dataDesignCowanRadio")$setActive(1)
				theWidget("dataDesignMoreyRadio")$setActive(0)
			}else if(wommbatAnalysis$Ktype==1) {
				theWidget("dataDesignPashlerRadio")$setActive(1)
				theWidget("dataDesignCowanRadio")$setActive(0)
				theWidget("dataDesignMoreyRadio")$setActive(0)
			}else if(wommbatAnalysis$Ktype==2){
				theWidget("dataDesignPashlerRadio")$setActive(0)
				theWidget("dataDesignCowanRadio")$setActive(0)
				theWidget("dataDesignMoreyRadio")$setActive(1)
			}
			
		}
	
		# Load columns
		gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("models"))
		if(!is.null(wommbatAnalysis$effs))
		{
			# Disable the column selection
			.womActiveColumnSelection(FALSE)
	
			# Setup the models tab initially
			.womInitializeModelsPageColumns()
			.womPopulateEffectsTable(wommbatAnalysis$effs)
	
			# Enable the models tab
			.womActiveModelTab(TRUE)
		}
		# set models
		if(length(wommbatAnalysis$Models)>0)
		{	
			gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("analysis"))
			.womActiveAnalysisTab(TRUE)
			for(i in 1:length(wommbatAnalysis$Models)){
				myName = names(wommbatAnalysis$Models)[i]
				.womAddModelToTreeview(myName)
				if(!is.null(wommbatAnalysis$Models[[myName]]$results))
				{
					gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("diagnostics"))
					.womActiveDiagnosticsTab(TRUE)
					.womActiveResultsTab(TRUE)
					.womSetResultsInTableview(myName)
				}
			}
		}
		
	
	
	}
	
}

.womLoadProjectFile <- function(projectFile)
{
	myNewEnv = new.env(parent=globalenv())
	loaded = try(load(projectFile,envir=myNewEnv), silent=T)
	if (inherits(loaded, "try-error")){
		gWidgets::gmessage(paste("The file '",projectFile,"' could not be loaded. The probable cause of this error is that the file is not an R save file.",
						sep=""), title="File error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	
	.womLoadProject(myNewEnv$wommbatAnalysis,projectFile)	
}


.womCreateNewColumnsOfInterestTreeModel <- function()
{
	model2 <- gtkTreeStoreNew("gchararray","gchararray","gboolean")
	treeview2 <- theWidget("dataTreeview2")
	
	gtkTreeViewSetModel(treeview2, model2)
}


.womCreateNewColumnSelectionTreeModel <- function()
{
	model1 <- gtkTreeStoreNew("gchararray","gboolean")
	treeview1 <- theWidget("dataTreeview1")
	
	gtkTreeViewSetModel(treeview1, model1)
}

.womCreateColumnsOfInterestColumns <- function()
{
	treeview2 <- theWidget("dataTreeview2")
	model2 <- gtkTreeViewGetModel(treeview2)
	


	#Create right columns tree view
	if(is.null(gtkTreeViewGetColumn(treeview2,0))){
		
		# column name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview2$insertColumnWithAttributes(-1, "Column", renderer, 
  								text = .womColumnsOfInterestTreeCols("columnname"),
								visible = .womColumnsOfInterestTreeCols("visible"))
								
		column <- treeview2$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# column type
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview2$insertColumnWithAttributes(-1, "Type", renderer, 
  								text = .womColumnsOfInterestTreeCols("type"),
								visible = .womColumnsOfInterestTreeCols("visible"))
								
		column <- treeview2$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}

}


.womSetDataForColumnSelection <- function(dataset,removeModels=TRUE,filename="")
{
	# Confirm that it is ok to reset the interface.
	if(length(wommbatAnalysis$Models)>0)
	{
		confirm = gWidgets::gconfirm("This will remove all existing models and analyses.", title="Reset interface?",toolkit=gWidgets::guiToolkit("RGtk2"))
		if(!confirm) return(0)
	}
	
	# check to make sure there are column names
	cols <- colnames(dataset)
	if(any(is.null(cols)))
	{
		gWidgets::gmessage("The data set could not be loaded, because there were no column names.", title="Data error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	
	# check to make sure that there are enough columns and rows.
	myDims = dim(dataset)
	if(myDims[1]<2 | myDims[2]<4)
	{
		gWidgets::gmessage("The data set could not be loaded, because too few rows (<2) or columns (<4).", title="Data error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)

	}
	
	# check to make sure all column names are unique. If not, error.
	if(any(table(cols))>1){
	gWidgets::gmessage("The data set could not be loaded, because some column names were duplicates.", title="Data error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	wommbatAnalysis$data <- dataset
	theWidget("dataFilenameEntry")$setText(filename)
	
	.womCreateNewColumnsOfInterestTreeModel()
	.womCreateNewColumnSelectionTreeModel()
	.womCreateColumnSelectionColumns()
	.womCreateColumnsOfInterestColumns()
	
	treeview1 <- theWidget("dataTreeview1")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 <- theWidget("dataTreeview2")
	model2 <- gtkTreeViewGetModel(treeview2)
	

	StateEnv$ColNamesIters=list()
	
	# Add dataset column names to treeview model
	for(col in cols){
		iter <- model1$append(NULL)$iter
		model1$set(iter,
	  		  .womColumnsTreeCols("columnname"), col,
			  .womColumnsTreeCols("visible"), TRUE)
		StateEnv$ColNamesIters[[col]] = iter
	}
	
	theWidget("dataRespEntry")$setText("")
	theWidget("dataChangeEntry")$setText("")
	theWidget("dataSetsizeEntry")$setText("")
	
	.womActiveColumnSelection(TRUE)
	.womActiveSaveTab(TRUE)
	.womClearInterfaceForDataLoad(removeModels)
	#StateEnv$win$present()	
}

.womClearInterfaceForDataLoad <- function(removeModels=TRUE)
{
	# Models
	if(removeModels) wommbatAnalysis$Models=NULL
	
	if(!is.null(theWidget("modelsEffectsTreeview")$getModel()))
		{
			theWidget("modelsEffectsTreeview")$getModel()$clear()
		}
	if(!is.null(theWidget("modelsOnKTreeview")$getModel()))
		{
			theWidget("modelsOnKTreeview")$getModel()$clear()
		}
	if(!is.null(theWidget("modelsOnZTreeview")$getModel()))
		{
			theWidget("modelsOnZTreeview")$getModel()$clear()
		}		
	if(!is.null(theWidget("modelsOnGTreeview")$getModel()))
		{
			theWidget("modelsOnGTreeview")$getModel()$clear()
		}
	if(!is.null(theWidget("modelsCovTreeview")$getModel()))
		{
			theWidget("modelsCovTreeview")$getModel()$clear()
		}
	if(!is.null(theWidget("modelsDefinedModelsTreeview")$getModel()))
		{
			theWidget("modelsDefinedModelsTreeview")$getModel()$clear()
		}		
	.womActiveModelTab(FALSE)
	theWidget("modelsNameEntry")$setText("Model")
	
	# Analysis
	if(!is.null(theWidget("analysisDefinedModelsTreeview")$getModel()))
		{
			theWidget("analysisDefinedModelsTreeview")$getModel()$clear()
		}		
	.womActiveAnalysisTab(FALSE)
	
	# Diagnostics
	
	if(!is.null(theWidget("diagnosticDefinedModelsTreeview")$getModel()))
		{
			theWidget("diagnosticDefinedModelsTreeview")$setModel(NULL)
		}		
	.womActiveDiagnosticsTab(FALSE)
	
	# Diagnostic plots
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
	# Results
	if(!is.null(theWidget("resultsDefinedModelsTreeview")$getModel()))
		{
			theWidget("resultsDefinedModelsTreeview")$setModel(NULL)	
		}		
	if(!is.null(theWidget("resultsParEstTreeview")$getModel()))
		{
			theWidget("resultsParEstTreeview")$setModel(NULL)
		}		
	
	.womActiveResultsTab(FALSE)
	
}

.womActiveColumnSelection<-function(status=TRUE)
{

	theWidget("dataColumnsFrame")$setSensitive(status)
	theWidget("dataBottomHbox")$setSensitive(status)

}

.womGetColumnNameSelection<-function()
{
	treeview = theWidget("dataTreeview1")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	vis <- model$get(iter, .womColumnsTreeCols("visible"))[[1]]
	return(list(iter=iter,visible=vis))
}


.womGetColumnOfInterestSelection<-function()
{
	treeview = theWidget("dataTreeview2")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	return(list(iter=iter))
}



.wom_clicked_saved_analysis <- function(button)
{
	#freezeGUI()
	#on.exit(thawGUI())
	fileChoose(action=".womLoadProjectFile", 
		type="open", 
		text="Select an Rdata file...", 
		filter = list("R data files" = list(patterns = c("*.Rdata","*.RData")),
						"All files" = list(patterns = c("*"))
					 )
		)
	
	
	#filename <- try(file.choose(), silent=T)
	#StateEnv$win$present()
	#if (inherits(filename, "try-error")) return()
	
	#.womLoadProjectFile(filename)
}



.wom_clicked_data_setsize <- function(button)
{
	treeview = theWidget("dataTreeview1")
	model <- gtkTreeViewGetModel(treeview)
	selection = .womGetColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model$get(selection$iter, .womColumnsTreeCols("columnname"))[[1]]

	
	# Check column to see if it is all integers
	colnum = match(selCol,colnames(wommbatAnalysis$data))
		if(!is.na(colnum)){
			if(!is.integer(wommbatAnalysis$data[,colnum]) | any(wommbatAnalysis$data[,colnum]<1))
			{
				# Inappropriate data in column
				gWidgets::gmessage(paste("Invalid set size column data:",selCol), title="Data error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
				return(0)
			}
		}else{
			gWidgets::gmessage("Could not find selected column in data set! This should not have happened.",
				title="Interface error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		}

		
	
	
	# If we had something there before, make it visible again
	oldText =  theWidget("dataSetsizeEntry")$getText()
	if(oldText!=""){
			model$set(StateEnv$ColNamesIters[[oldText]], .womColumnsTreeCols("visible"), TRUE)
	}
	
	# Set the column name in the entry box
	theWidget("dataSetsizeEntry")$setText(selCol)
	
	# Set the selected iter to be invisible
	model$set(selection$iter, .womColumnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
	
}

.wom_clicked_data_response <- function(button)
{
	treeview = theWidget("dataTreeview1")
	model <- gtkTreeViewGetModel(treeview)
	selection = .womGetColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model$get(selection$iter, .womColumnsTreeCols("columnname"))[[1]]

	
	# Check column to see if it is all zeros and ones
	colnum = match(selCol,colnames(wommbatAnalysis$data))
		if(!is.na(colnum)){
			if(!identical(sort(as.integer(unique(wommbatAnalysis$data[,colnum]))),as.integer(c(0,1))))
			{
				# Inappropriate data in column
				gWidgets::gmessage(paste("Invalid response column data:", selCol), title="Data error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
				return(0)
			}
		}else{
			gWidgets::gmessage("Could not find selected column in data set! This should not have happened.",
				title="Interface error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(0)
		}

	
	# If we had something there before, make it visible again
	oldText =  theWidget("dataRespEntry")$getText()
	if(oldText!=""){
			model$set(StateEnv$ColNamesIters[[oldText]], .womColumnsTreeCols("visible"), TRUE)
	}

	# Set the column name in the entry box
	theWidget("dataRespEntry")$setText(selCol)
	
	# Set the selected iter to be invisible
	model$set(selection$iter, .womColumnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
	
}

.wom_clicked_data_change <- function(button)
{
	treeview = theWidget("dataTreeview1")
	model <- gtkTreeViewGetModel(treeview)
	selection = .womGetColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model$get(selection$iter, .womColumnsTreeCols("columnname"))[[1]]

	
	# Check column to see if it is all zeros and ones
	colnum = match(selCol,colnames(wommbatAnalysis$data))
		if(!is.na(colnum)){
			if(!identical(sort(as.integer(unique(wommbatAnalysis$data[,colnum]))),as.integer(c(0,1))))
			{
				# Inappropriate data in column
				gWidgets::gmessage(paste("Invalid change column data:",selCol), title="Data error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
				return(0)
			}
		}else{
			gWidgets::gmessage("Could not find selected column in data set! This should not have happened.",
				title="Interface error",icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(0)
		}

	
	# If we had something there before, make it visible again
	oldText =  theWidget("dataChangeEntry")$getText()
	if(oldText!=""){
			model$set(StateEnv$ColNamesIters[[oldText]], .womColumnsTreeCols("visible"), TRUE)
	}

	# Set the column name in the entry box
	theWidget("dataChangeEntry")$setText(selCol)
		
	# Set the selected iter to be invisible
	model$set(selection$iter, .womColumnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
	
}

.wom_clicked_data_categ <- function(button)
{
	.womAddColumnOfInterest("categorical")
}

.wom_clicked_data_contin <- function(button)
{
	.womAddColumnOfInterest("continuous")
}

.womAddColumnOfInterest <-function(type)
{
	treeview1 = theWidget("dataTreeview1")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 = theWidget("dataTreeview2")
	model2 <- gtkTreeViewGetModel(treeview2)
	
	selection = .womGetColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model1$get(selection$iter, .womColumnsTreeCols("columnname"))[[1]]
	
	# Add the column name in treeview2
	iter <- model2$append(NULL)$iter
	model2$set(iter, 
			  .womColumnsOfInterestTreeCols("columnname"), selCol,
			  .womColumnsOfInterestTreeCols("type"), type,
			  .womColumnsOfInterestTreeCols("visible"), TRUE)
		
	
		
	# Set the selected iter to be invisible
	model1$set(selection$iter, .womColumnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
}

.wom_clicked_data_remove <- function(button)
{
	treeview1 = theWidget("dataTreeview1")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 = theWidget("dataTreeview2")
	model2 <- gtkTreeViewGetModel(treeview2)
	
	selection = .womGetColumnOfInterestSelection()
	
	# If nothing is selected, do nothing
	if(is.null(selection$iter)){
		return()
	}
	
	# Get the name of the column, make visible again
	selCol <- model2$get(selection$iter, .womColumnsOfInterestTreeCols("columnname"))[[1]]
	iter1 <- StateEnv$ColNamesIters[[selCol]]
	model1$set(iter1, .womColumnsTreeCols("visible"), TRUE)
	
	# Remove from treeview2
	model2$remove(selection$iter)
	
	# Set the selected iter to be invisible
	StateEnv$win$present()
}

.wom_clicked_unlock_data <-function(button)
{
	if(!is.null(wommbatAnalysis$data))
	{
		.womSetDataForColumnSelection(wommbatAnalysis$data)
	}
}


.wom_clicked_get_data_frame <- function(button)
{
	freezeGUI()
	on.exit(thawGUI())
	.womSelectDataFrame(envir = .GlobalEnv)
}




.wom_clicked_data_apply <- function(button)
{
	treeview2 = theWidget("dataTreeview2")
	model2 <- gtkTreeViewGetModel(treeview2)

	# Check to make sure we defined at least one column of interest
	nRows <- gtkTreeModelIterNChildren(model2, NULL)	
	if(nRows==0)
	{
		gWidgets::gmessage(paste("You must define at least one column of interest."), title="Data error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	
	# Check to make sure we've defined necessary columns
	respText <- wommbatAnalysis$response <- theWidget("dataRespEntry")$getText()
	changeText <- wommbatAnalysis$change <- theWidget("dataChangeEntry")$getText()
	setsizeText <- wommbatAnalysis$setsize <- theWidget("dataSetsizeEntry")$getText()
	
	
	# Get model type
	if(theWidget("dataDesignMoreyRadio")$getActive()){
		wommbatAnalysis$Ktype=as.integer(2)
	}
	if(theWidget("dataDesignPashlerRadio")$getActive()){
		wommbatAnalysis$Ktype=as.integer(1)
	}
	if(theWidget("dataDesignCowanRadio")$getActive()){
		wommbatAnalysis$Ktype=as.integer(0)
	}
	# Get model type
	
	
	if(respText=="" | changeText=="" | setsizeText=="")
	{
		gWidgets::gmessage(paste("You must define the response, change, and setsize columns."), title="Data error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)
	}
	
	# Disable the column selection
	.womActiveColumnSelection(FALSE)
	
	# Setup the models tab initially
	.womInitializeModelsPageColumns()
	.womInitialModelSetup()
	
	# Enable the models tab
	.womActiveModelTab(TRUE)
	
	# Turn the page
	gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("models"))

	
	StateEnv$win$present()
}

