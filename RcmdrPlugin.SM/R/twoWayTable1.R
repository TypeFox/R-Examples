
twoWayTable1 <- function(){ # dialog memory 2011-06-27 J. Fox
	Library("abind")
	Library("RColorBrewer")
	Library("vcd")
	defaults <- list(initial.row=NULL, initial.column=NULL, 
			initial.percents="row", 
			initial.subset=gettextRcmdr("<all valid cases>"))
	dialog.values <- getDialog("twoWayTable1", defaults)
	initializeDialog(title=gettextRcmdr("Graphiques (C*C)"))
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
		putDialog("twoWayTable1", list(
						initial.row=row, 
						initial.column=column, 
						initial.percents=percents,initial.subset=initial.subset
				))
		if (length(row) == 0 || length(column) == 0){
			errorCondition(recall=twoWayTable1, message=gettextRcmdr("You must select two variables."))
			return()
		}
		if (row == column) {
			errorCondition(recall=twoWayTable1, message=gettextRcmdr("Row and column variables are the same."))
			return()
		}
		closeDialog()
		

		command0 <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
				subset, ")", sep="")
		

		
		if (percents == "row") {
command<-paste("mosaicplot(",command0,",main=\"\",col=brewer.pal(max(length(dimnames(",command0,")[[2]]),3),\"Set1\"))",sep="")
doItAndPrint(command)
}
		if (percents == "column"){
command<-paste("barplot(t(sweep(",command0,",1,apply(",command0,",1,sum),\"/\")),beside=TRUE,col=brewer.pal(max(length(dimnames(",command0,")[[2]]),3),\"Set1\")[1:(length(dimnames(",command0,")[[2]]))],xlab=names(dimnames(",command0,"))[1],legend.text=TRUE,ylab=\"pourcentage\")",sep="")
doItAndPrint(command)}

		if (percents == "total"){
command<-paste("assoc(",command0,",shade=TRUE)",sep="")
doItAndPrint(command)}
		
		

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="xtabs", reset="twoWayTable1")
	radioButtons(name="percents",
			buttons=c("rowPercents", "columnPercents", "totalPercents"),
			values=c("row", "column", "total"), initialValue=dialog.values$initial.percents,
			labels=gettextRcmdr(c("mosaique", "en barres", "association")), title=gettextRcmdr("Type de graphique"))
	
	tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
	tkgrid(variablesFrame, sticky="w")
	tkgrid(percentsFrame, sticky="w")
	
	
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}