
twoWayTable3 <- function(){ # dialog memory 2011-06-27 J. Fox
	Library("abind")

	defaults <- list(initial.row=NULL, initial.column=NULL, 
			initial.percents="row", 
			initial.subset=gettextRcmdr("<all valid cases>"))
	dialog.values <- getDialog("twoWayTable3", defaults)
	initializeDialog(title=gettextRcmdr("Significativit\xE9 (C*C)"))
	variablesFrame <- tkframe(top)
	.factors <- Factors()
	rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
			initialSelection=varPosn(dialog.values$initial.row, "factor"))
	columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
			initialSelection=varPosn(dialog.values$initial.column, "factor"))
	subsetBox(subset.expression=dialog.values$initial.subset)
	onOK <- function(){
		row <- getSelection(rowBox)
		column <- getSelection(columnBox)
		
		
		
		
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
				else paste(", subset=", subset, sep="")
		putDialog("twoWayTable3", list(
						initial.row=row, 
						initial.column=column, 
						initial.subset=initial.subset
				))
		if (length(row) == 0 || length(column) == 0){
			errorCondition(recall=twoWayTable3, message=gettextRcmdr("You must select two variables."))
			return()
		}
		if (row == column) {
			errorCondition(recall=twoWayTable3, message=gettextRcmdr("Row and column variables are the same."))
			return()
		}
		closeDialog()
		

		command <- paste("signifCC(xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
				subset, "))", sep="")
		
		doItAndPrint(command)
		
		
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="xtabs", reset="twoWayTable3")
		
	tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
	tkgrid(variablesFrame, sticky="w")
	
	
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}