# based on rcmdr code enterTable

enterTableMHTest <- function(){
  Library("abind")
  env <- environment()
  defaults <- list(initial.tab=0)
  dialog.values <- getDialog("enterTableMHTest", defaults)
  initializeDialog(title=gettextRcmdr("Enter Two-Way Table for Marginal Homogeneity test"), use.tabs=TRUE, tabs=c("tableTab", "statisticsTab"))
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
  initial.test <- if (is.null(initial.table)) "default" else attr(initial.table, "test")
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
    putDialog("enterTableMHTest", list(initial.tab=tab))
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
      errorCondition(recall=enterTableMHTest, message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
      return()
    }
    if (length(unique(row.names)) != nrows){
      errorCondition(recall=enterTableMHTest, message=gettextRcmdr("Row names are not unique."))
      return()
    }
    if (length(unique(col.names)) != ncols){
      errorCondition(recall=enterTableMHTest, message=gettextRcmdr("Column names are not unique."))
      return()
    }
    percents <- as.character(tclvalue(percentsVariable))
    test <- as.character(tclvalue(testVariable)) #
    
    closeDialog()
    
    
#     command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
#                      ", byrow=TRUE)", sep="")
#     doItAndPrint(paste(".Table <- ", command, sep=""))    
    command <- paste("local({\n  .Table <- matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols, ", byrow=TRUE)
                      \n  cat('\\nFrequency table:\\n')\n  print(.Table)", sep="")

    if (percents == "row") 
      command <- paste(command, '\n  cat("\\nRow percentages:\\n")\n  print(rowPercents(.Table))',
                       sep="")
    else if (percents == "column") 
      command <-  paste(command, '\n  cat("\\nColumn percentages:\\n")\n  print(colPercents(.Table))',
                        sep="")
    else if (percents == "total") 
      command <- paste(command, '\n  cat("\\nTotal percentages:\\n")\n  print(totPercents(.Table))',
                       sep="")
    
    if (test == "default") {
      strDistribution = ""
    } else {
      strDistribution = paste(", distribution='", test,"' ", sep="")
    } 
    
    command <- paste(command, "\n  .Test <- mh_test(as.table(.Table)", strDistribution,")", sep="")
    #command.2 <- paste(command.2, "\n  .Test <- mh_test(as.table(.Table)", strDistribution,")", sep="")
    command <- paste(command, "\n  print(.Test)", sep="")
    
    command <- paste(command, "\n})", sep="")
    doItAndPrint(command)
    
    
    
    
    

#     command <- paste("c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), ")", sep="")
#     justDoIt(paste("rownames(.Table) <- ", command, sep=""))
#     logger(paste("rownames(.Table) <- ", command, sep=""))
#     command <- paste("c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), ")", sep="")
#     justDoIt(paste("colnames(.Table) <- ", command, sep=""))
#     logger(paste("colnames(.Table) <- ", command, sep=""))
#     doItAndPrint(".Table  # Counts")
#     if (percents == "row") doItAndPrint("rowPercents(.Table) # Row Percentages")
#     if (percents == "column") doItAndPrint("colPercents(.Table) # Column Percentages")
#     if (percents == "total") doItAndPrint("totPercents(.Table) # Percentage of Total")
#          if (test == "default") {
#            strDistribution = ""
#          } else {
#            strDistribution = paste(", distribution='", test,"' ", sep="")
#          } 
#         
#         command <- paste(".Test <- mh_test(as.table(.Table)", strDistribution,")", sep="")
#         doItAndPrint(command)
#         command <- paste(".Test", sep="")
#         doItAndPrint(command)
#     
#         if (getRcmdr("retain.selections")){
#           attr(.Table, "percentages") <- percents
#           #attr(.Table, "test") <- "default"
#           attr(.Table, "test") <- test
#           putRcmdr("savedTableCoin", .Table)
#         }
#     
#         logger("remove(.Test)")
#         remove(.Test, envir=.GlobalEnv)
#     
#         logger("remove(.Table)")
#         remove(.Table, envir=.GlobalEnv)
    if (getRcmdr("retain.selections")){
      attr(.Table, "percentages") <- percents
      attr(.Table, "tests") <- c(chisq, chisqComp, expected, fisher)
      putRcmdr("savedTable", .Table)
    }
    logger("remove(.Table)")
    remove(.Table, envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="mh_test", reset="resetenterTableMHTest", apply = "enterTableMHTest")
  radioButtons(statisticsTab, name="percents", buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"), values=c("row", "column", "total", "none"),
               initialValue=initial.percentages, labels=gettextRcmdr(c("Row percentages", "Column percentages",  "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
#   checkBoxes(statisticsTab, frame="testsFrame", boxes=c("chisq", "chisqComponents", "expFreq", "fisher"), initialValues=initial.tests,
#              labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
#                                    "Print expected frequencies", "Fisher's exact test")), title=gettextRcmdr("Hypothesis Test"))
  tkgrid(labelRcmdr(sliderFrame, text=gettextRcmdr("Number of Rows:")), rowsSlider, rowsShow, sticky="we", padx = 6,  pady = 6)
  tkgrid(labelRcmdr(sliderFrame, text=gettextRcmdr("Number of Columns:")), colsSlider, colsShow, sticky="we", padx = 6,  pady = 6)
  tkgrid(sliderFrame, sticky="w")
  tkgrid(labelRcmdr(tableTab, text=gettextRcmdr("Enter counts:"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="we", padx = 6,  pady = 6)
  tkgrid(percentsFrame, sticky="we", padx = 6,  pady = 6)
  radioButtons(statisticsTab, name="test", buttons=c("default", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test")) 
  tkgrid(testFrame, sticky="w")
#   tkgrid(testsFrame, sticky="we", padx = 6, pady = 6)
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("tableTab", "statisticsTab"), tab.names=c("Table", "Statistics"))
}
resetenterTableMHTest <- function(){
  putRcmdr("savedTable", NULL)
  enterTableMHTest()
}



fncCoinMHTest <- function(){
  Library("abind")
  defaults <- list(initial.row=NULL, initial.column=NULL, initial.test="default", 
                   initial.percents="none", initial.subset=gettextRcmdr("<all valid cases>"), initial.tab=0)
  dialog.values <- getDialog("fncCoinMHTest", defaults)
  initializeDialog(title=gettextRcmdr("Marginal Homogeneity Test for Two-Way Table"), use.tabs=TRUE)
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
    
    test <- as.character(tclvalue(testVariable)) #
    
    percents <- as.character(tclvalue(percentsVariable))
    initial.subset <- subset <- tclvalue(subsetVariable)
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
    else paste(", subset=", subset, sep="")
    putDialog("fncCoinMHTest", list(
      initial.row=row, 
      initial.column=column, 
      initial.test=test, 
      initial.percents=percents, initial.subset=initial.subset,
      initial.tab=tab))
      if (length(row) == 0 || length(column) == 0){
        errorCondition(recall=fncCoinMHTest, message=gettextRcmdr("You must select two variables."))
        return()
      }
      if (row == column) {
        errorCondition(recall=fncCoinMHTest, message=gettextRcmdr("The variables are the same."))
        return()
      }
    
    command <- paste("length(levels(as.factor(", ActiveDataSet(), "$",row, "))) == length(levels(as.factor(", ActiveDataSet(), "$", column, ")))", sep="")
    command.2 <- paste("local({\nputRcmdr('.bolEqualNbLevels', ", command, ")\n})")
    justDoIt(command.2)
    .bolEqualNbLevels2 <- getRcmdr(".bolEqualNbLevels")
      if (.bolEqualNbLevels2 == FALSE) {
        errorCondition(recall=fncCoinMHTest, message=gettextRcmdr("The factors have different number of levels."))
        return()
      }
      
      command <- paste("levels(as.factor(", ActiveDataSet(), "$",row, ")) == levels(as.factor(", ActiveDataSet(), "$", column, "))", sep="")
      command.2 <- paste("local({\nputRcmdr('.bolEqualLevels', ", command, ")\n})")
      justDoIt(command.2)
      .bolEqualLevels2 <- getRcmdr(".bolEqualLevels")    
      if (.bolEqualLevels2 == FALSE) {
        errorCondition(recall=fncCoinMHTest, message=gettextRcmdr("The factors have different levels."))
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
    
    if (test == "default") {
      strDistribution = ""
    } else {
      strDistribution = paste(", distribution='", test,"' ", sep="")
    } 
    
    command <- paste(command, "\n  .Test <- mh_test(as.table(.Table)", strDistribution,")", sep="")
    command.2 <- paste(command.2, "\n  .Test <- mh_test(as.table(.Table)", strDistribution,")", sep="")
    command <- paste(command, "\n  print(.Test)", sep="")
    
    command <- paste(command, "\n})", sep="")
    doItAndPrint(command)
    
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="mh_test", reset="fncCoinMHTest", apply="fncCoinMHTest")
  radioButtons(optionsTab, name="percents",
               buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
               values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
               labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), 
               title=gettextRcmdr("Compute Percentages"))
  tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
  
  
  tkgrid(variablesFrame, sticky="w")
  tkgrid(percentsFrame, sticky="w")
  radioButtons(optionsTab, name="test", buttons=c("default", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test")) 
  tkgrid(testFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tab.names=c("Data", "Statistics"))
}