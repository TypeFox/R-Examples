frequencyDistribution2 <- function () {
	defaults <- list (initial.x = NULL)
	dialog.values <- getDialog ("frequencyDistribution2", defaults)
	initializeDialog(title = gettextRcmdr("Frequency Distributions"))
	xBox <- variableListBox(top, Factors(), selectmode = "multiple", 
			title = gettextRcmdr("Variables (pick one or more)"),
			initialSelection = varPosn (dialog.values$initial.x, "factor"))
	optionsFrame <- tkframe(top)
	

	onOK <- function() {
		x <- getSelection(xBox)
		if (length(x) == 0) {
			errorCondition(recall = frequencyDistribution2, message = gettextRcmdr("You must select a variable."))
			return()
		}
		
		putDialog ("frequencyDistribution2", list (initial.x = x))
		
		closeDialog()
		.activeDataSet <- ActiveDataSet()
		for (variable in x) {
			command <- paste("table(", .activeDataSet, "$", variable, 
					")# effectifs pour ",variable, sep = "")
			
			
			doItAndPrint(command)


			command <- paste("percent(table(", .activeDataSet, "$", variable, 
					"))# pourcentages pour ", 
							variable, sep = "")
			
			doItAndPrint(command)
		}
		
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "table", reset = "frequencyDistribution2")
	tkgrid(getFrame(xBox), sticky = "nw")
	
	tkgrid(optionsFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 3, columns = 2)
}

