barCode <- function () {
	Library("barcode")
	defaults <- list (initial.variable = NULL)
	dialog.values <- getDialog ("barCode", defaults)
	initializeDialog(title = gettextRcmdr("Graphe en code barre"))
	variableBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.variable, "numeric"))



	onOK <- function() {
		variable <- getSelection(variableBox)


		
putDialog ("barCode", list (initial.variable = variable))

		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall = barCode, message = gettextRcmdr("You must select a variable"))
			return()
		}
		.activeDataSet <- ActiveDataSet()


		command <- (paste("barcode(", .activeDataSet, "$", 
							variable, ",xlab=\"",variable,"\")", sep = ""))



		
		logger(command)
		justDoIt(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "barcode", reset = "barCode")
	tkgrid(getFrame(variableBox), sticky = "nw")


	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 3, columns = 1)
}
