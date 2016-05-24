
bbbivariate0<- function () {
	
	defaults <- list(initial.group = NULL, initial.response = NULL)
	dialog.values <- getDialog("bbbivariate0", defaults)
	initializeDialog(title = gettextRcmdr(paste("Graphiques crois", "\U00E9", "s en saucisson",sep = "")))
	groupBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variable explicative"), 
			initialSelection = varPosn(dialog.values$initial.group, "all"))
	responseBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"),selectmode = "multiple",
			initialSelection = varPosn(dialog.values$initial.response, "all"))
	
	

	
onOK <- function() {
	
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
putDialog ("bbbivariate0", list (initial.group = group, initial.response=response))

		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = bbbivariate0, message = gettextRcmdr("You must select a groups factor."))
			return()
		}
		if (length(response) < 2) {
			errorCondition(recall = bbbivariate0, message = gettextRcmdr("You must select at least two response variables."))
			return()
		}

		.activeDataSet <- ActiveDataSet()
listvar <- paste(response, collapse = "\",\"")

command <- paste("bivariate0(",.activeDataSet,"[,c(\"", 
            listvar, "\")],",.activeDataSet,"$",group,",xlab=\"",group,"\")",sep = "")

logger(command)
justDoIt(command)
		activateMenus()
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "bivariate0", model = TRUE, reset = "bbbivariate0")
	
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	
tkgrid(buttonsFrame, sticky = "w")



	dialogSuffix(rows = 2, columns = 2)
}
