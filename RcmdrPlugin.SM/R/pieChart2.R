pieChart2 <- function () {
	Library("colorspace")
	defaults <- list (initial.variable = NULL,initial.mean="0")
	dialog.values <- getDialog ("pieChart2", defaults)
	initializeDialog(title = gettextRcmdr("Pie Chart"))
	variableBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.variable, "factor"))

checkBoxes(frame="checkBoxFrame", boxes=c("mean"), 
               initialValues=c(dialog.values$initial.mean), 
               labels=gettextRcmdr(c("Pareto (rangement en fonction des effectifs)")))

	onOK <- function() {
		variable <- getSelection(variableBox)

meanVar <- tclvalue(meanVariable)
		
putDialog ("pieChart", list (initial.variable = variable, initial.mean=meanVar))

		closeDialog()
		if (length(variable) == 0) {
			errorCondition(recall = pieChart2, message = gettextRcmdr("You must select a variable"))
			return()
		}
		.activeDataSet <- ActiveDataSet()

if (meanVar == 0){
		command <- (paste("pie(table(", .activeDataSet, "$", 
							variable, "), labels=levels(", .activeDataSet, "$", 
							variable, "), main=\"", variable, "\", col=rainbow_hcl(length(levels(", .activeDataSet, "$", variable,")),start=30,end=-300))", sep = ""))
}
else
{
command <- (paste("pie(rev(sort(table(", .activeDataSet, "$", 
							variable, "))), main=\"", variable, "\", col=heat_hcl(length(levels(", .activeDataSet, "$", variable,")),c=c(80,30),l=c(30,90),power=c(1/5,2)),clockwise=TRUE)", sep = ""))
}


		
		logger(command)
		justDoIt(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "pie", reset = "pieChart2")
	tkgrid(getFrame(variableBox), sticky = "nw")

tkgrid(checkBoxFrame, sticky="w")

	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix(rows = 3, columns = 1)
}
