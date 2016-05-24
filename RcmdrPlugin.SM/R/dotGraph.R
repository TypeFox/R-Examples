dotGraph <- function () {
# change	
defaults <- list (initial.variable = NULL,initial.mean="0")
	
dialog.values <- getDialog ("dotGraph", defaults)
	initializeDialog(title = gettextRcmdr("Dot Graph"))
	variableBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.variable, "factor"))
# change
checkBoxes(frame="checkBoxFrame", boxes=c("mean"), 
               initialValues=c(dialog.values$initial.mean), 
               labels=gettextRcmdr(c("Rangement en fonction des effectifs (Pareto)")))

	onOK <- function() {
		variable <- getSelection(variableBox)
# change
meanVar <- tclvalue(meanVariable)

#change
		putDialog ("dotGraph", list (initial.variable = variable, initial.mean=meanVar))
		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall = dotGraph, message = gettextRcmdr("You must select a variable"))
			return()
		}

#change
if (meanVar == 0){
		command <- paste("dotchartTable(table(", ActiveDataSet(), "$", 
				variable, "),pch=16, ylab=\"", variable, "\", xlab=\"Effectif\")", 
				sep = "")
}
else
{
command <- paste("dotchartTable(sort(table(", ActiveDataSet(), "$", 
				variable, ")), pch=16,ylab=\"", variable, "\", xlab=\"Effectif\")", 
				sep = "")
}

		logger(command)
		justDoIt(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "dotchart", reset = "dotGraph")
	tkgrid(getFrame(variableBox), sticky = "nw")

#change
tkgrid(checkBoxFrame, sticky="w")

	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 2, columns = 1)
}





