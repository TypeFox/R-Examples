# this stores application (non-project) state information
if (!exists("StateEnv", environment(), inherits=FALSE)) {
	StateEnv <- new.env()
}
# stores all project information
if (!exists("wommbatAnalysis", environment(), inherits=FALSE)) {
	wommbatAnalysis <- new.env()
}

wommbatGUI <- function(project = NULL, projectFile= NULL,CSVfile = NULL, dataFrame = NULL, devel = FALSE) {

	options("guiToolkit"="RGtk2")

	StateEnv$update <- list()
	wommbatAnalysis$Models <- list()
	
	StateEnv$echo.to.console <- TRUE
	StateEnv$echo.to.log <- TRUE
	StateEnv$Graphics <- list()
	
	#StateEnv$GUI <- gladeXMLNew(getpackagefile("wommbatGlade.glade"),
	#	root="window1")

	StateEnv$GUI <- gtkBuilder()
	filename <- getpackagefile("wommbatMain.ui")
	res <- StateEnv$GUI$addFromFile(filename)
	if (!is.null(res$error))
		stop("ERROR: ", res$error$message)
		
	StateEnv$win <- theWidget("window1")
	StateEnv$win$setTitle("WoMMBAT")
	
	# connect the callbacks (event handlers)
	#gladeXMLSignalAutoconnect(StateEnv$GUI)
	StateEnv$GUI$connectSignals(NULL)
	StateEnv$handlers = list()
	
	# Connect scrollbar signal
	scrollBar <- theWidget("hscrollbar1")
	adjustment <- gtkAdjustmentNew(value = 1, lower = 1, upper = 1, step.incr=1)
	scrollBar$setAdjustment(adjustment)
	gtkAdjustmentSetValue(scrollBar$getAdjustment(), 1)
	StateEnv$handlers$diagnosticScrollBar1 <- gSignalConnect(theWidget("hscrollbar1"), "value-changed", .scrolled_diagnostics_scrollbar)
	
	
	resultsSelection <- theWidget("resultsParEstTreeview")$getSelection()		
	gtkTreeSelectionSetMode(resultsSelection, "GTK_SELECTION_MULTIPLE")
	gSignalConnect(resultsSelection, "changed", .resultsSelectionChanged)

	.womSetupDiagnosticsComboBox()	
	.loadTooltips()
	.populateTextviews()
	
	.womSetInitialSensitivity()
	
	# load files/dataframes
	if(!is.null(dataFrame)){
		if(is.data.frame(dataFrame)){
			theWidget("dataFilenameEntry")$setText("<Loaded from data.frame>")
			.womSetDataForColumnSelection(dataFrame)
		}
	}else if(!is.null(CSVfile))
	{
		.womOpenCSVFile(CSVfile)
	}else if(!is.null(project))
	{
		.womLoadProject(project)
	}else if(!is.null(projectFile))
	{
		.womLoadProjectFile(projectFile)
	}
	
	StateEnv$win$present()
	if(!devel){
		theWidget("dataDesignMoreyRadio")$hide()
	}
	return(invisible(NULL))
}


wommbatNoGUI <- function(project=NULL, projectFile= NULL, settings)
{
	stop('This function is not yet implemented.')
	return(invisible(NULL))
}

