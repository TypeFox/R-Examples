
enConstruction <- function () {
	
	defaults <- list(initial.variables = NULL) 
	dialog.values <- getDialog("Enconstruction", defaults)
	
	initializeDialog(title = gettextRcmdr("Attention Travaux !"))
	
	onOK <- function() {
		
		closeDialog()
		

			command <- "print(\"un peu de patience\")"
			#logger(command)
			doItAndPrint(command)
		
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "dotchart", reset = "plotLikert")
	
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 6, columns = 2)
}

