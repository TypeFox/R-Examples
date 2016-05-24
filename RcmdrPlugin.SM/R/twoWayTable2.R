
twoWayTable2 <- function(){ # dialog memory 2011-06-27 J. Fox
	Library("abind")
	defaults <- list(initial.row=NULL, initial.column=NULL, 
			initial.percents="none", 
			initial.subset=gettextRcmdr("<all valid cases>"))
	dialog.values <- getDialog("twoWayTable2", defaults)
	initializeDialog(title=gettextRcmdr("Two-Way Table"))
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
		percents <- as.character(tclvalue(percentsVariable))
		
		
		
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
				else paste(", subset=", subset, sep="")
		putDialog("twoWayTable2", list(
						initial.row=row, 
						initial.column=column, 
						initial.percents=percents,initial.subset=initial.subset
				))
		if (length(row) == 0 || length(column) == 0){
			errorCondition(recall=twoWayTable2, message=gettextRcmdr("You must select two variables."))
			return()
		}
		if (row == column) {
			errorCondition(recall=twoWayTable2, message=gettextRcmdr("Row and column variables are the same."))
			return()
		}
		closeDialog()
		command0 <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
				subset, ")", sep="")
		
		doItAndPrint(command0)
		if (percents == "row") {
command<-paste("rowPercents(",command0,") # Pourcentages en lignes",sep="")
doItAndPrint(command)
}

		if (percents == "column"){
command<-paste("colPercents(",command0,") # Pourcentages en colonnes",sep="")
doItAndPrint(command)
} 


		if (percents == "total") {
command<-paste("totPercents(",command0,") # Pourcentages du Total",sep="")
doItAndPrint(command)
} 


		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="xtabs", reset="twoWayTable2")
	radioButtons(name="percents",
			buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
			values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
			labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
	
	tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
	tkgrid(variablesFrame, sticky="w")
	tkgrid(percentsFrame, sticky="w")
	
	
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}