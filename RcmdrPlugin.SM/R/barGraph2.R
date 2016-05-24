barGraph2 <- function () {
Library("colorspace")	
defaults <- list (initial.variable = NULL,initial.mean="0")

	dialog.values <- getDialog ("barGraph2", defaults)
	initializeDialog(title = gettextRcmdr("Bar Graph"))
	variableBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.variable, "factor"))

checkBoxes(frame="checkBoxFrame", boxes=c("mean"), 
               initialValues=c(dialog.values$initial.mean), 
               labels=gettextRcmdr(c("Pareto (rangement en fonction des effectifs)")))

	onOK <- function() {
		variable <- getSelection(variableBox)

meanVar <- tclvalue(meanVariable)

		putDialog ("barGraph", list (initial.variable = variable,initial.mean=meanVar))
		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall = barGraph2, message = gettextRcmdr("You must select a variable"))
			return()
		}

if (meanVar == 0){
		command <- paste("barplot(table(", ActiveDataSet(), "$", 
				variable, "), xlab=\"", variable, "\", ylab=\"Effectif\",col=\"pink\")", 
				sep = "")
}
else
{
command <- paste("barplot(rev(sort(table(", ActiveDataSet(), "$", 
				variable, "))), xlab=\"", variable, "\", ylab=\"Effectif\",col=heat_hcl(length(levels(", .activeDataSet, "$", variable,")),c=c(80,30),l=c(30,90),power=c(1/5,2)))", 
				sep = "")
}



		
		logger(command)
		justDoIt(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "barplot", reset = "barGraph2")
	tkgrid(getFrame(variableBox), sticky = "nw")

tkgrid(checkBoxFrame, sticky="w")

	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 2, columns = 1)
}



