mmultivariate0<- function () {
Library("ade4")

	
	defaults <- list(initial.listx = NULL, initial.listy = NULL)
	dialog.values <- getDialog("mmultivariate0", defaults)
	initializeDialog(title = gettextRcmdr(paste("Analyse multivari", "\U00E9", "e factorielle", 
            sep = "")))
	xBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables explicatives"),selectmode = "multiple",
			initialSelection = varPosn(dialog.values$initial.listx, "all"))
	yBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables supplementaires"),selectmode = "multiple",
			initialSelection = varPosn(dialog.values$initial.listy, "all"))
	
	

	
onOK <- function() {
	
		listx <- getSelection(xBox)
		listy <- getSelection(yBox)
putDialog ("mmultivariate0", list (initial.listx = listx, initial.listy=listy))

		closeDialog()
		if (length(listx) < 2) {
			errorCondition(recall = mmultivariate0, message = gettextRcmdr("You must select at least 2 explanatory variables."))
			return()
		}

		.activeDataSet <- ActiveDataSet()
listvarx <- paste(listx, collapse = "\",\"")

if (length(listy)==0){
command <- paste("multivariate0(",.activeDataSet,"[,c(\"", 
            listvarx, "\")])",sep = "")
}
else
{
listvary <- paste(listy, collapse = "\",\"")
command <- paste("multivariate0(",.activeDataSet,"[,c(\"", 
            listvarx, "\")],",.activeDataSet,"[,c(\"", 
            listvary,"\")])",sep = "")
}
logger(command)
justDoIt(command)

		activateMenus()
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "multivariate0", model = TRUE, reset= "mmultivariate0")
	
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	
tkgrid(buttonsFrame, sticky = "w")



	dialogSuffix(rows = 2, columns = 2)
}



