# based on code from Rcmdr by J. Fox twoWayTable

fncCoinTwoWayTable <- function(){
  Library("abind")
  defaults <- list(initial.row=NULL, initial.column=NULL, 
                   initial.percents="none", initial.chisq=1, initial.chisqComp=0, initial.expected=0, 
                   initial.fisher=0, initial.adjPval=0, initial.lblTest=0, initial.subset=gettextRcmdr("<all valid cases>"), initial.tab=0)
  dialog.values <- getDialog("fncCoinTwoWayTable", defaults)
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
    adjPvalues <- tclvalue(adjPvalVariable)
    linearbylinear <- tclvalue(lblTestVariable)
    fisher <- tclvalue(fisherTestVariable)
    initial.subset <- subset <- tclvalue(subsetVariable)
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
    else paste(", subset=", subset, sep="")
    putDialog("fncCoinTwoWayTable", list(
      initial.row=row, 
      initial.column=column, 
      initial.percents=percents, initial.chisq=chisq, initial.chisqComp=chisqComp, 
      initial.expected=expected, initial.adjPval=adjPvalues, initial.lblTest=linearbylinear, initial.fisher=fisher, initial.subset=initial.subset,
      initial.tab=tab))
    if (length(row) == 0 || length(column) == 0){
      errorCondition(recall=fncCoinTwoWayTable, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (row == column) {
      errorCondition(recall=fncCoinTwoWayTable, message=gettextRcmdr("Row and column variables are the same."))
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
      if (adjPvalues == 1) {
        command <- paste(command, '\n  cat("\\np-values adjusted for multiple testing by a single-step max-T multiple testing approach:\\n")\n  print(pvalue(independence_test(', column, ' ~ ', row, ', teststat = "max", data=', ActiveDataSet(),
                     subset, '), method = "single-step"))', sep="")
        #doItAndPrint(paste("pvalue(independence_test(", column, " ~ ", row, ", teststat = 'max', data=", ActiveDataSet(),
        #             subset, "), method = 'single-step') # p-values adjusted for multiple testing by a single-step max-T multiple testing approach", sep=''))
      }
    }
    if (linearbylinear == 1) {
      command <- paste(command, "\n  print(lbl_test(.Table))") 
      #command.2 <- paste("local({\n  .warn <- options(warn=-1)\nprint(lbl_test(.Table))\n  options(.warn)\n})")
      #justDoIt(command.2)
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
  OKCancelHelp(helpSubject="chisq_test", reset="fncCoinTwoWayTable", apply="fncCoinTwoWayTable")
  radioButtons(optionsTab, name="percents",
               buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
               values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
               labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), 
               title=gettextRcmdr("Compute Percentages"))
  checkBoxes(optionsTab, frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "adjPval", "lblTest", "fisherTest"), 
             initialValues=c(dialog.values$initial.chisq, dialog.values$initial.chisqComp, 
                             dialog.values$initial.expected, dialog.values$initial.adjPval, dialog.values$initial.lblTest, dialog.values$initial.fisher),
             labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
                                   "Print expected frequencies", "Print adjusted p-values for standardized contingency table","Linear-by-linear association test", "Fisher's exact test")))
  tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(percentsFrame, sticky="w")
  tkgrid(labelRcmdr(optionsTab, text=gettextRcmdr("Hypothesis Tests"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
  tkgrid(testsFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tab.names=c("Data", "Statistics"))
}

# fncCoinTwoWayTable <- function(){
#   Library("abind")
#   initializeDialog(title=gettextRcmdr("Two-Way Table"))
#   variablesFrame <- tkframe(top)
#   .factors <- Factors()
#   rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable\n(select one)"))
#   columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable\n(select one)"))
#   subsetBox()
#   onOK <- function(){
#     row <- getSelection(rowBox)
#     column <- getSelection(columnBox)
#     if (length(row) == 0 || length(column) == 0){
#       errorCondition(recall=fncCoinTwoWayTable, message=gettextRcmdr("You must select two variables."))
#       return()
#     }
#     if (row == column) {
#       errorCondition(recall=fncCoinTwoWayTable, message=gettextRcmdr("Row and column variables are the same."))
#       return()
#     }
#     percents <- as.character(tclvalue(percentsVariable))
#     chisq <- tclvalue(chisqTestVariable)
#     chisqComp <- tclvalue(chisqComponentsVariable)
#     expected <- tclvalue(expFreqVariable)
#     adjPvalues <- tclvalue(adjPvalVariable)
#     linearbylinear <- tclvalue(lblTestVariable)
#     fisher <- tclvalue(fisherTestVariable)
#     subset <- tclvalue(subsetVariable)
#     subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
#     else paste(", subset=", subset, sep="")
#     closeDialog()
#     command <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
#                      subset, ")", sep="")
#     logger(paste(".Table <- ", command, sep=""))
#     assign(".Table", justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(".Table")
#     if (percents == "row") doItAndPrint("rowPercents(.Table) # Row Percentages")
#     if (percents == "column") doItAndPrint("colPercents(.Table) # Column Percentages")
#     if (percents == "total") doItAndPrint("totPercents(.Table) # Percentage of Total")
#     if (chisq == 1) {
#       command <- "chisq.test(.Table, correct=FALSE)" # using old test to get info for expected counts
#       logger(paste(".Test <- ", command, sep=""))
#       assign(".Test", justDoIt(command), envir=.GlobalEnv)
#       #doItAndPrint(".Test") not showing clasical chi square test
#       doItAndPrint("chisq_test(.Table)") #showing goin chi square test
#       if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
#       warnText <- NULL
#       if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
#                                                                    gettextRcmdr("expected frequencies are less than 1"))
#       if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
#                                                                    gettextRcmdr(" expected frequencies are less than 5"), sep="")
#       if (!is.null(warnText)) Message(message=warnText,
#                                       type="warning")
#       if (chisqComp == 1) {
#         command <- "round(.Test$residuals^2, 2) # Chi-square Components"
#         doItAndPrint(command)
#       }
#       if (adjPvalues == 1) {
#         doItAndPrint(paste("pvalue(independence_test(", column, " ~ ", row, ", teststat = 'max', data=", ActiveDataSet(),
#                      subset, "), method = 'single-step') # p-values adjusted for multiple testing by a single-step max-T multiple testing approach", sep=''))
#       }
#       logger("remove(.Test)")
#       remove(.Test, envir=.GlobalEnv)
#     }
#     if (linearbylinear == 1) doItAndPrint("lbl_test(.Table)")
#     if (fisher == 1) doItAndPrint("fisher.test(.Table)")
#     logger("remove(.Table)")
#     remove(.Table, envir=.GlobalEnv)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="chisq_test")
#   radioButtons(name="percents",
#                buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
#                values=c("row", "column", "total", "none"), initialValue="none",
#                labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
#   checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "adjPval", "lblTest", "fisherTest"), initialValues=c("1", "0", "0", "0", "0", "0"),
#              labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
#                                    "Print expected frequencies", "Print adjusted p-values for standardized contingency table","Linear-by-linear association test", "Fisher's exact test")))
#   tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
#   tkgrid(variablesFrame, sticky="w")
#   tkgrid(percentsFrame, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
#   tkgrid(testsFrame, sticky="w")
#   tkgrid(subsetFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=7, columns=1)
# }