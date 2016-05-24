bbbisummarize0<- function () {
Library("ade4")
	
	defaults <- list(initial.group = NULL, initial.response = NULL,initial.mean="0")
	dialog.values <- getDialog("bbbisummarize0", defaults)
	initializeDialog(title = gettextRcmdr(paste("Statistiques crois", "\U00E9", "es en saucisson",sep = "")))
	groupBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variable explicative"), 
			initialSelection = varPosn(dialog.values$initial.group, "all"))
	responseBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"),selectmode = "multiple",
			initialSelection = varPosn(dialog.values$initial.response, "all"))
checkBoxes(frame="checkBoxFrame", boxes=c("mean"), 
               initialValues=c(dialog.values$initial.mean), 
               labels=gettextRcmdr(c("Rangement des variables en fonction de la p-value")))	
	

	
onOK <- function() {
	
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
		meanVar <- tclvalue(meanVariable)

putDialog ("bbbisummarize0", list (initial.group = group, initial.response=response,initial.mean=meanVar))

		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = bbbisummarize0, message = gettextRcmdr("You must select a explanatory variable."))
			return()
		}
		if (length(response) <2) {
			errorCondition(recall = bbbisummarize0, message = gettextRcmdr("You must select at least select two response variables."))
			return()
		}

		.activeDataSet <- ActiveDataSet()
listvar <- paste(response, collapse = "\",\"")

if (meanVar == 0){
command <- paste("bisummarize0(",.activeDataSet,"[,c(\"", 
            listvar, "\")],",.activeDataSet,"$",group,",tri=FALSE)",sep = "")
}
else{
command <- paste("bisummarize0(",.activeDataSet,"[,c(\"", 
            listvar, "\")],",.activeDataSet,"$",group,",tri=TRUE)",sep = "")
}


doItAndPrint(command)


		activateMenus()
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "bisummarize0", model = TRUE, reset = "bbbisummarize0")
	
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	tkgrid(checkBoxFrame, sticky="w")

tkgrid(buttonsFrame, sticky = "w")



	dialogSuffix(rows = 2, columns = 2)
}


