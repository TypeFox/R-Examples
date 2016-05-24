
plotLikert <- function () {
	
	defaults <- list(initial.variables = NULL) 
	dialog.values <- getDialog("plotLikert", defaults)
	
	initializeDialog(title = gettextRcmdr("Graphe des moyennes"))
	variablesBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Select variables (three or more)"), 
			selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "numeric"))
	onOK <- function() {
		variables <- getSelection(variablesBox)
		closeDialog()
		if (length(variables) < 3) {
			errorCondition(recall =plotLikert, message = gettextRcmdr("Fewer than 3 variable selected."))
			return()
		}
				putDialog("plotLikert", list(initial.variables = variables))
		.activeDataSet <- ActiveDataSet() 


		listvar<-paste(variables,collapse="\",\"")

			command <- paste("dotchart(sort(sapply(",.activeDataSet,"[c(\"",listvar,"\")],mean,na.rm=TRUE)),xlab=\"Moyennes\",pch=16)", sep = "")
			logger(command)
			justDoIt(command)
		
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot)
	OKCancelHelp(helpSubject = "dotchart", reset = "plotLikert")
	tkgrid(getFrame(variablesBox), sticky = "nw")
	
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 6, columns = 2)
}

