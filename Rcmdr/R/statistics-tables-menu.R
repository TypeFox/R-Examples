# Statistics Menu dialogs

# last modified 2014-08-01 by J. Fox

# Tables menu

twoWayTable <- function(){
    Library("abind")
    defaults <- list(initial.row=NULL, initial.column=NULL, 
        initial.percents="none", initial.chisq=1, initial.chisqComp=0, initial.expected=0, 
        initial.fisher=0, initial.subset=gettextRcmdr("<all valid cases>"), initial.tab=0)
    dialog.values <- getDialog("twoWayTable", defaults)
    initializeDialog(title=gettextRcmdr("Two-Way Table"), use.tabs=TRUE)
    variablesFrame <- tkframe(dataTab)
    .factors <- Factors()
    rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
        initialSelection=varPosn(dialog.values$initial.row, "factor"))
    columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
        initialSelection=varPosn(dialog.values$initial.column, "factor"))
    subsetBox(dataTab, subset.expression=dialog.values$initial.subset)
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
        initial.subset <- subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
        else paste(", subset=", subset, sep="")
        putDialog("twoWayTable", list(
            initial.row=row, 
            initial.column=column, 
            initial.percents=percents, initial.chisq=chisq, initial.chisqComp=chisqComp, 
            initial.expected=expected, initial.fisher=fisher, initial.subset=initial.subset,
            initial.tab=tab))
        if (length(row) == 0 || length(column) == 0){
            errorCondition(recall=twoWayTable, message=gettextRcmdr("You must select two variables."))
            return()
        }
        if (row == column) {
            errorCondition(recall=twoWayTable, message=gettextRcmdr("Row and column variables are the same."))
            return()
        }
        closeDialog()
        command <- paste("local({\n  .Table <- xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
            subset, ')\n  cat("\\nFrequency table:\\n")\n  print(.Table)', sep="")
        command.2 <- paste("local({\n  .warn <- options(warn=-1)\n  .Table <- xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
                           subset, ")", sep="")
        if (percents == "row") 
          command <- paste(command, '\n  cat("\\nRow percentages:\\n")\n  print(rowPercents(.Table))',
                           sep="")
        else if (percents == "column") 
          command <-  paste(command, '\n  cat("\\nColumn percentages:\\n")\n  print(colPercents(.Table))',
                            sep="")
        else if (percents == "total") 
          command <- paste(command, '\n  cat("\\nTotal percentages:\\n")\n  print(totPercents(.Table))',
                sep="")
        if (chisq == 1) {
            command <- paste(command, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
            command.2 <- paste(command.2, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
            command <- paste(command, "\n  print(.Test)", sep="")
            if (expected == 1)
                command <- paste(command, '\n  cat("\\nExpected counts:\\n")\n  print(.Test$expected)', sep="")
            if (chisqComp == 1) {
                command <- paste(command, '\n  cat("\\nChi-square components:\\n")\n  print(round(.Test$residuals^2, 2))', sep="")
            }
        }
        if (fisher == 1) command <- paste(command, "\n  print(fisher.test(.Table))")
        command <- paste(command, "\n})", sep="")
        doItAndPrint(command)
        if (chisq == 1){
          command.2 <- paste(command.2, "\nputRcmdr('.expected.counts', .Test$expected)\n  options(.warn)\n})")
          justDoIt(command.2)
          warnText <- NULL
          expected <- getRcmdr(".expected.counts")
          if (0 < (nlt1 <- sum(expected < 1))) warnText <- paste(nlt1,
                                                                       gettextRcmdr("expected frequencies are less than 1"))
          if (0 < (nlt5 <- sum(expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                                                                       gettextRcmdr(" expected frequencies are less than 5"), sep="")
          if (!is.null(warnText)) Message(message=warnText,
                                          type="warning")
        }
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="xtabs", reset="twoWayTable", apply="twoWayTable")
    radioButtons(optionsTab, name="percents",
        buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
        values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), 
        title=gettextRcmdr("Compute Percentages"))
    checkBoxes(optionsTab, frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), 
        initialValues=c(dialog.values$initial.chisq, dialog.values$initial.chisqComp, 
            dialog.values$initial.expected, dialog.values$initial.fisher),
        labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
            "Print expected frequencies", "Fisher's exact test")))
    tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(optionsTab, text=gettextRcmdr("Hypothesis Tests"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tab.names=c("Data", "Statistics"))
}

multiWayTable <- function (){
	Library("abind")
	defaults <- list (initial.row = NULL, initial.column = NULL, initial.control = NULL, 
			initial.percents = "none", initial.subset=gettextRcmdr("<all valid cases>"))
	dialog.values <- getDialog ("multiWayTable", defaults)
	initializeDialog(title = gettextRcmdr("Multi-Way Table"))
	variablesFrame <- tkframe(top)
	.factors <- Factors()
	rowBox <- variableListBox(variablesFrame, .factors, title = gettextRcmdr("Row variable (pick one)"),
			initialSelection = varPosn (dialog.values$initial.row, "factor"))
	columnBox <- variableListBox(variablesFrame, .factors, title = gettextRcmdr("Column variable (pick one)"),
			initialSelection = varPosn (dialog.values$initial.column, "factor"))
	controlBox <- variableListBox(variablesFrame, .factors, selectmode = "multiple", 
			title = gettextRcmdr("Control variable(s) (pick one or more)"), 
			initialSelection = varPosn (dialog.values$initial.control, "factor"))
	subsetBox(subset.expression = dialog.values$initial.subset)
	onOK <- function() {
		row <- getSelection(rowBox)
		column <- getSelection(columnBox)
		controls <- getSelection(controlBox)
		if (length(row) == 0 || length(column) == 0 || length(controls) == 
				0) {
			errorCondition(recall = multiWayTable, message = gettextRcmdr("You must select row, column, and control variables"))
			return()
		}
		if ((row == column) || is.element(row, controls) || is.element(column, 
				controls)) {
			errorCondition(recall = multiWayTable, message = gettextRcmdr("Row, column, and control variables must be different."))
			return()
		}
		percents <- as.character(tclvalue(percentsVariable))
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) 
					""
				else paste(", subset=", subset, sep = "")
		putDialog ("multiWayTable", list (initial.row = row, initial.column = column, initial.control = controls, initial.percents = percents, initial.subset=initial.subset))
		closeDialog()
		command <- paste("local({\n  .Table <- xtabs(~", row, "+", column, "+", paste(controls, 
						collapse = "+"), ", data=", ActiveDataSet(), subset, 
				')\n  cat("\\nFrequency table:\\n")\n  print(.Table)', sep = "")
		if (percents == "row") 
			command <- paste(command, '\n  cat("\\nRow percentages:\\n")\n  print(rowPercents(.Table))', sep="")
		if (percents == "column") 
		  command <- paste(command, '\n  cat("\\nColumn percentages:\\n")\n  print(colPercents(.Table))', sep="")
	  command <- paste(command, "\n})")
    doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "xtabs", reset = "multiWayTable", apply = "multiWayTable")
	radioButtons(name = "percents", buttons = c("rowPercents", 
					"columnPercents", "nonePercents"), values = c("row", 
					"column", "none"),  labels = gettextRcmdr(c("Row percentages", 
							"Column percentages", "No percentages")), title = gettextRcmdr("Compute Percentages"),
			initialValue = dialog.values$initial.percents)
	tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text = "    "), 
			getFrame(columnBox), labelRcmdr(variablesFrame, text = "    "), 
			getFrame(controlBox), sticky = "nw")
	tkgrid(variablesFrame, sticky = "w")
	tkgrid(percentsFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "ew")
	dialogSuffix()
}

enterTable <- function(){
  Library("abind")
  env <- environment()
  defaults <- list(initial.tab=0)
  dialog.values <- getDialog("enterTable", defaults)
  initializeDialog(title=gettextRcmdr("Enter Two-Way Table"), use.tabs=TRUE, tabs=c("tableTab", "statisticsTab"))
  assign(".tableFrame", tkframe(tableTab), envir=env)  
  setUpTable <- function(...){
    tkdestroy(get(".tableFrame", envir=env))
    assign(".tableFrame", tkframe(tableTab), envir=env)
    nrows <- as.numeric(tclvalue(rowsValue))
    ncols <- as.numeric(tclvalue(colsValue))
    make.col.names <- "labelRcmdr(.tableFrame, text='')"
    for (j in 1:ncols) {
      col.varname <- paste(".colname.", j, sep="")
      assign(col.varname, if (is.null(initial.table) || j > length(colnames)) tclVar(j) else tclVar(colnames[j]), envir=env)
      make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                              col.varname, ")", sep="")
    }
    eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
    for (i in 1:nrows){
      varname <- paste(".tab.", i, ".1", sep="")
      assign(varname, if (is.null(initial.table) || i > length(rownames)) tclVar("") else tclVar(initial.table[i, 1]) , envir=env)
      row.varname <- paste(".rowname.", i, sep="")
      assign(row.varname, if (is.null(initial.table) || i > length(rownames)) tclVar(i) else tclVar(rownames[i]), envir=env)
      make.row <- paste("ttkentry(.tableFrame, width='5', textvariable=",
                        row.varname, ")", sep="")
      make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                        varname, ")", sep="")
      for (j in 2:ncols){
        varname <- paste(".tab.", i, ".", j, sep="")
        assign(varname, if (is.null(initial.table) || i > length(rownames) || j > length(colnames)) 
          tclVar("") else tclVar(initial.table[i, j]), envir=env)
        make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                          varname, ")", sep="")
      }
      eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
    }
    tkgrid(get(".tableFrame", envir=env), sticky="ew", padx = 6)
  }
  initial.table <- getRcmdr("savedTable")
  initial.percentages <- if (is.null(initial.table)) "none" else attr(initial.table, "percentages")
  initial.tests <- if (is.null(initial.table)) c("1", "0", "0", "0") else attr(initial.table, "tests")
  if (is.null(initial.table)){
    rowsValue <- tclVar("2")
    colsValue <- tclVar("2")
  }
  else {
    rowsValue <- tclVar(nrow(initial.table))
    colsValue <- tclVar(ncol(initial.table))
    rownames <- rownames(initial.table)
    colnames <- colnames(initial.table)
  }
  sliderFrame <- tkframe(tableTab)
  rowsSlider <- tkscale(sliderFrame, from=2, to=10, showvalue=FALSE, variable=rowsValue,
                        resolution=1, orient="horizontal", command=setUpTable)
  rowsShow <- labelRcmdr(sliderFrame, textvariable=rowsValue, width=2, justify="right")
  colsSlider <- tkscale(sliderFrame, from=2, to=10, showvalue=FALSE, variable=colsValue,
                        resolution=1, orient="horizontal", command=setUpTable)
  colsShow <- labelRcmdr(sliderFrame, textvariable=colsValue, width=2, justify="right")
  onOK <- function(){
    tab <- if (as.character(tkselect(notebook)) == tableTab$ID) 0 else 1
    putDialog("enterTable", list(initial.tab=tab))
    nrows <- as.numeric(tclvalue(rowsValue))
    ncols <- as.numeric(tclvalue(colsValue))
    cell <- 0
    counts <- rep(NA, nrows*ncols)
    row.names <- rep("", nrows)
    col.names <- rep("", ncols)
    for (i in 1:nrows) row.names[i] <-
      eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
    for (j in 1:ncols) col.names[j] <-
      eval(parse(text=paste("tclvalue(", paste(".colname.", j, sep=""),")", sep="")))
    for (i in 1:nrows){
      for (j in 1:ncols){
        cell <- cell+1
        varname <- paste(".tab.", i, ".", j, sep="")
        counts[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
      }
    }
    counts <- na.omit(counts)
    if (length(counts) != nrows*ncols){
      errorCondition(recall=enterTable, message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
      return()
    }
    if (length(unique(row.names)) != nrows){
      errorCondition(recall=enterTable, message=gettextRcmdr("Row names are not unique."))
      return()
    }
    if (length(unique(col.names)) != ncols){
      errorCondition(recall=enterTable, message=gettextRcmdr("Column names are not unique."))
      return()
    }
    percents <- as.character(tclvalue(percentsVariable))
    chisq <- tclvalue(chisqVariable)
    chisqComp <- tclvalue(chisqComponentsVariable)
    expected <- tclvalue(expFreqVariable)
    fisher <- tclvalue(fisherVariable)
    closeDialog()
    command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
                     ", byrow=TRUE)", sep="")
    doItAndPrint(paste(".Table <- ", command, sep=""))
    command <- paste("c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), ")", sep="")
    justDoIt(paste("rownames(.Table) <- ", command, sep=""))
    logger(paste("rownames(.Table) <- ", command, sep=""))
    command <- paste("c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), ")", sep="")
    justDoIt(paste("colnames(.Table) <- ", command, sep=""))
    logger(paste("colnames(.Table) <- ", command, sep=""))
    doItAndPrint(".Table  # Counts")
    if (percents == "row") doItAndPrint("rowPercents(.Table) # Row Percentages")
    if (percents == "column") doItAndPrint("colPercents(.Table) # Column Percentages")
    if (percents == "total") doItAndPrint("totPercents(.Table) # Percentage of Total")
    if (chisq == 1) {
      command <- "chisq.test(.Table, correct=FALSE)"
      doItAndPrint(paste(".Test <- ", command, sep=""))
      doItAndPrint(".Test")
      if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
      warnText <- NULL
      if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                                                                   gettextRcmdr("expected frequencies are less than 1"))
      if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                                                                   gettextRcmdr(" expected frequencies are less than 5"), sep="")
      if (!is.null(warnText)) Message(message=warnText,
                                      type="warning")
      if (chisqComp == 1) {
        command <- "round(.Test$residuals^2, 2) # Chi-square Components"
        doItAndPrint(command)
      }
      logger("remove(.Test)")
      remove(.Test, envir=.GlobalEnv)
    }
    if (fisher == 1) doItAndPrint("fisher.test(.Table)")
    if (getRcmdr("retain.selections")){
      attr(.Table, "percentages") <- percents
      attr(.Table, "tests") <- c(chisq, chisqComp, expected, fisher)
      putRcmdr("savedTable", .Table)
    }
    logger("remove(.Table)")
    remove(.Table, envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="chisq.test", reset="resetEnterTable", apply = "enterTable")
  radioButtons(statisticsTab, name="percents", buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"), values=c("row", "column", "total", "none"),
               initialValue=initial.percentages, labels=gettextRcmdr(c("Row percentages", "Column percentages",  "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
  checkBoxes(statisticsTab, frame="testsFrame", boxes=c("chisq", "chisqComponents", "expFreq", "fisher"), initialValues=initial.tests,
             labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
                                   "Print expected frequencies", "Fisher's exact test")), title=gettextRcmdr("Hypothesis Test"))
  tkgrid(labelRcmdr(sliderFrame, text=gettextRcmdr("Number of Rows:")), rowsSlider, rowsShow, sticky="we", padx = 6,  pady = 6)
  tkgrid(labelRcmdr(sliderFrame, text=gettextRcmdr("Number of Columns:")), colsSlider, colsShow, sticky="we", padx = 6,  pady = 6)
  tkgrid(sliderFrame, sticky="w")
  tkgrid(labelRcmdr(tableTab, text=gettextRcmdr("Enter counts:"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="we", padx = 6,  pady = 6)
  tkgrid(percentsFrame, sticky="we", padx = 6,  pady = 6)
  tkgrid(testsFrame, sticky="we", padx = 6, pady = 6)
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("tableTab", "statisticsTab"), tab.names=c("Table", "Statistics"))
}
resetEnterTable <- function(){
    putRcmdr("savedTable", NULL)
    enterTable()
}
