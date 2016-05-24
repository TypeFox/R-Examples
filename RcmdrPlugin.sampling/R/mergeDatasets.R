# New merge function
# Same as the old but have taken out forced by=rownames line

mergeDataSets <- function(){
	dataSets <- listDataSets()
	.activeDataSet <- ActiveDataSet()
	initializeDialog(title=gettextRcmdr("Merge Data Sets"))
	dsname <- tclVar("MergedDataset")
	dsnameFrame <- tkframe(top)
	entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
	dataSet1Box <- variableListBox(top, dataSets, title=gettextRcmdr("First Data Set (pick one)"),
		initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
	dataSet2Box <- variableListBox(top, dataSets, title=gettextRcmdr("Second Data Set (pick one)"))
	commonVar <- tclVar("0")
	commonFrame <- tkframe(top)
	commonButton <- tkcheckbutton(commonFrame, variable=commonVar)	
	radioButtons(top, "direction", buttons=c("rows", "columns"), 
		labels=gettextRcmdr(c("Merge rows", "Merge columns")), title=gettextRcmdr("Direction of Merge"))
	onOK <- function(){
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == "") {
			errorCondition(recall=mergeDataSets,
				message=gettextRcmdr("You must enter the name of a data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)) {
			errorCondition(recall=mergeDataSets,
				message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
				closeDialog()
				mergeDataSets()
				return()
			}
		}
		name1 <- getSelection(dataSet1Box)
		name2 <- getSelection(dataSet2Box)
		if (length(name1) == 0){
			errorCondition(recall=mergeDataSets,
				message=gettextRcmdr("You must select a data set."))
			return()
		}
		if (length(name2) == 0){
			errorCondition(recall=mergeDataSets,
				message=gettextRcmdr("You must select a data set."))
			return()
		}
		if (name1 == name2){
			errorCondition(recall=mergeDataSets,
				message=gettextRcmdr("You cannot merge a data set with itself."))
			return()
		}
		common <- if (tclvalue(commonVar) == "1") TRUE else FALSE
		direction <- tclvalue(directionVariable)
		if (direction == "rows"){
			command <- paste(dsnameValue, " <- mergeRows(", name1, ", ", name2,
				", common.only=", common, ")", sep="")
			doItAndPrint(command)	
		}
		else {
			command <- paste(dsnameValue, " <- merge(", name1, ", ", name2,
				", all=", !common, sep="")
			doItAndPrint(command)
		}
		activeDataSet(dsnameValue)
		closeDialog()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="mergeRows")
	tkgrid(labelRcmdr(dsnameFrame, text=gettextRcmdr("Name for merged data set:  ")), entryDsname)
	tkgrid(dsnameFrame, sticky="w", columnspan=2)
	tkgrid(getFrame(dataSet1Box), getFrame(dataSet2Box), sticky="nw")
	tkgrid(labelRcmdr(commonFrame, text=gettextRcmdr("Merge only common\nrows or columns")), 
		commonButton, sticky="nw")
	tkgrid(directionFrame, commonFrame, sticky="sw")
	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=5, columns=2)
}