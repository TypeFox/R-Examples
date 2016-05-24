ppyramide0 <- function () {
Library("Hmisc")
	defaults <- list(initial.group = NULL, initial.response = NULL)
	dialog.values <- getDialog("ppyramide0", defaults)
	initializeDialog(title = gettextRcmdr("Pyramide des ages"))
	groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable Sexe"), 
			initialSelection = varPosn(dialog.values$initial.group, "factor"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable Age"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	
	

	
onOK <- function() {
	
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
putDialog ("ppyramide0", list (initial.group = group, initial.response=response))

		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = ppyramide0, message = gettextRcmdr("You must select a groups factor."))
			return()
		}
		if (length(response) == 0) {
			errorCondition(recall = ppyramide0, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		command <- paste("pyramide0(",.activeDataSet,"$",response,",",.activeDataSet,"$",group,")",sep="")
		logger(command)
		justDoIt(command)
		activateMenus()
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "pyramide0", model = TRUE, reset = "ppyramide0")
	
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	
tkgrid(buttonsFrame, sticky = "w")



	dialogSuffix(rows = 4, columns = 2)
}