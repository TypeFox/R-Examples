# based on code from Rcmdr by J. Fox multiWayTable

# fncCoinCochraneMHTest <- function(){
#   Library("abind")
#   initializeDialog(title="Cochran-Mantel-Haenzsel test")
#   variablesFrame <- tkframe(top)
#   .factors <- Factors()
#   rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable\n(select one)"))
#   columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable\n(select one)"))
#   controlBox <- variableListBox(variablesFrame, .factors, selectmode="multiple",
#                                 title=gettextRcmdr("Control variable(s)\n(select one or more)"))
#   subsetBox()
#   onOK <- function(){
#     row <- getSelection(rowBox)
#     column <- getSelection(columnBox)
#     controls <- getSelection(controlBox)
#     if (length(row) == 0 || length(column) == 0 || length(controls) == 0) {
#       errorCondition(recall=fncCoinCochraneMHTest , message=gettextRcmdr("You must select row, column, and control variables"))
#       return()
#     }
#     if ((row == column) || is.element(row, controls) || is.element(column, controls)) {
#       errorCondition(recall=fncCoinCochraneMHTest , message=gettextRcmdr("Row, column, and control variables must be different."))
#       return()
#     }
#     percents <- as.character(tclvalue(percentsVariable))
#     subset <- tclvalue(subsetVariable)
#     subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
#     else paste(", subset=", subset, sep="")
#     closeDialog()
#     command <- paste("xtabs(~", row, "+", column, "+", paste(controls, collapse="+"),
#                      ", data=", ActiveDataSet(), subset, ")", sep="")
#     logger(paste(".Table <- ", command, sep=""))
#     assign(".Table", justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(".Table")
#     if (percents == "row") doItAndPrint("rowPercents(.Table) # Row Percentages")
#     if (percents == "column") doItAndPrint("colPercents(.Table) # Column Percentages")
#     doItAndPrint(".Table")
#     doItAndPrint("cmh_test(.Table)")
#     remove(.Table, envir=.GlobalEnv)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="cmh_test")
#   radioButtons(name="percents", buttons=c("rowPercents", "columnPercents", "nonePercents"), values=c("row", "column", "none"),
#                initialValue="none", labels=gettextRcmdr(c("Row percentages", "Column percentages", "No percentages")), title=gettextRcmdr("Compute Percentages"))
#   tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), labelRcmdr(variablesFrame, text="    "),
#          getFrame(controlBox), sticky="nw")
#   tkgrid(variablesFrame, sticky="w")
#   tkgrid(percentsFrame, sticky="w")
#   tkgrid(subsetFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=4, columns=1)
# }


fncCoinCochraneMHTest <- function(){
  Library("abind")
  defaults <- list(initial.row=NULL, initial.column=NULL, initial.control=NULL, 
                   initial.percents="none", initial.subset=gettextRcmdr("<all valid cases>"), initial.tab=0)
  dialog.values <- getDialog("fncCoinCochraneMHTest", defaults)
  initializeDialog(title=gettextRcmdr("Cochran-Mantel-Haenzsel test"), use.tabs=TRUE)
  variablesFrame <- tkframe(dataTab)
  .factors <- Factors()
  rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
                            initialSelection=varPosn(dialog.values$initial.row, "factor"))
  columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
                               initialSelection=varPosn(dialog.values$initial.column, "factor"))
  controlBox <- variableListBox(variablesFrame, .factors, selectmode="multiple",
                                   title=gettextRcmdr("Control variable(s)\n(select one or more)"),
                                initialSelection=varPosn(dialog.values$initial.control, "factor"))
  subsetBox(dataTab, subset.expression=dialog.values$initial.subset)
  onOK <- function(){
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    row <- getSelection(rowBox)
    column <- getSelection(columnBox)
    
    controls <- getSelection(controlBox)
    
    percents <- as.character(tclvalue(percentsVariable))
    initial.subset <- subset <- tclvalue(subsetVariable)
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
    else paste(", subset=", subset, sep="")
    putDialog("fncCoinCochraneMHTest", list(
      initial.row=row, 
      initial.column=column, 
      initial.control=controls, 
      initial.percents=percents, initial.subset=initial.subset,
      initial.tab=tab))
    if (length(row) == 0 || length(column) == 0 || length(controls) == 0) {
      errorCondition(recall=fncCoinCochraneMHTest , message=gettextRcmdr("You must select row, column, and control variables"))
      return()
    }
    if ((row == column) || is.element(row, controls) || is.element(column, controls)) {
      errorCondition(recall=fncCoinCochraneMHTest , message=gettextRcmdr("Row, column, and control variables must be different."))
      return()
    }
    closeDialog()
    
    command <- paste("local({\n  .Table <- xtabs(~", row, "+", column, "+", paste(controls, collapse="+"), ", data=", ActiveDataSet(),
                     subset, ')\n  cat("\\nFrequency table:\\n")\n  print(.Table)', sep="")
    command.2 <- paste("local({\n  .warn <- options(warn=-1)\n  .Table <- xtabs(~", row, "+", column, "+", paste(controls, collapse="+"), ", data=", ActiveDataSet(),
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

      command <- paste(command, "\n  .Test <- cmh_test(.Table)", sep="")
      command.2 <- paste(command.2, "\n  .Test <- cmh_test(.Table)", sep="")
      command <- paste(command, "\n  print(.Test)", sep="")

    command <- paste(command, "\n})", sep="")
    doItAndPrint(command)

    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cmh_test", reset="fncCoinCochraneMHTest", apply="fncCoinCochraneMHTest")
  radioButtons(optionsTab, name="percents",
               buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
               values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
               labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), 
               title=gettextRcmdr("Compute Percentages"))
  tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), labelRcmdr(variablesFrame, text="    "),
        getFrame(controlBox), sticky="nw")
  
  tkgrid(variablesFrame, sticky="w")
  tkgrid(percentsFrame, sticky="w")
  #tkgrid(labelRcmdr(optionsTab, text=gettextRcmdr("Hypothesis Tests"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
  tkgrid(subsetFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tab.names=c("Data", "Statistics"))
}