## dialog memory added 21020823 rmh
## based on Rcmdr/Rcmdr_1.8-4.tar.gz!Rcmdr/R/statistics-tables-menu.R by John Fox

## Tables menu

twoWayTable.HH <- function() {
   Library("abind")
   defaults <- list(initial.row=NULL, initial.column=NULL,
                    initial.rowPct=0, initial.colPct=0, initial.totPct=0,
                    initial.chisq=1, initial.chisqComp=1, initial.chiComp=0,
                    initial.expected=0,
                    initial.fisher=0, initial.subset=gettextRcmdr("<all valid cases>"))
   dialog.values <- getDialog("twoWayTable.HH", defaults)
   initializeDialog(title=gettextRcmdr("Two-Way Table"))
   variablesFrame <- tkframe(top)
   .factors <- Factors()
   rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
                             initialSelection=varPosn(dialog.values$initial.row, "factor"))
   columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
                                initialSelection=varPosn(dialog.values$initial.column, "factor"))
   subsetBox(subset.expression=dialog.values$initial.subset)
   subsetBox(subset.expression=dialog.values$initial.subset)
    onOK <- function() {
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
        rowPct <- tclvalue(rowPercentsVariable)
        colPct <- tclvalue(colPercentsVariable)
        totPct <- tclvalue(totPercentsVariable)
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        chiComp <- tclvalue(chiCompVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)

        initial.subset <- subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
        else paste(", subset=", subset, sep="")
        putDialog("twoWayTable.HH",
                  list(
                    initial.row=row,
                    initial.column=column,
                    initial.rowPct=rowPct, initial.colPct=colPct, initial.totPct=totPct,
                    initial.chisq=chisq, initial.chisqComp=chisqComp, initial.chiComp=chiComp,
                    initial.expected=expected, initial.fisher=fisher, initial.subset=initial.subset
                    ))
        if (length(row) == 0 || length(column) == 0){
          errorCondition(recall=twoWayTable.HH,
                         message=gettextRcmdr("You must select two variables."))
          return()
        }
        if (row == column) {
          errorCondition(recall=twoWayTable.HH,
                         message=gettextRcmdr("Row and column variables are the same."))
          return()
        }
        closeDialog()

        command <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
            subset, ")", sep="")
        justDoIt(paste(".Table <- ", command, sep=""))
        doItAndPrint(".Table")

        if (rowPct == 1) doItAndPrint("rowPercents(.Table) # Row Percentages")
        if (colPct == 1) doItAndPrint("colPercents(.Table) # Column Percentages")
        if (totPct == 1) doItAndPrint("totPercents(.Table) # Percentage of Total")

        if (chisq == 1) {
            command <- "chisq.test(.Table, correct=FALSE)"
            justDoIt(paste(".Test <- ", command, sep=""))
            doItAndPrint(".Test")
            if (chisqComp == 1) doItAndPrint("round(.Test$residuals^2, 2) # Chi-square Components")
            if (chiComp == 1) doItAndPrint("round(.Test$residuals, 2) # Chi Components (residuals)")
            if (expected == 1) doItAndPrint("round(.Test$expected, 2) # Expected Counts")
##            if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                gettextRcmdr("expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                gettextRcmdr(" expected frequencies are less than 5"), sep="")
            if (!is.null(warnText)) Message(message=warnText,
                type="warning")
            }
        if (fisher == 1) doItAndPrint("fisher.test(.Table)")
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="xtabs", reset="twoWayTable.HH")

    checkBoxes(frame="percentsFrame",
               boxes=c(
                 "rowPercents",
                 "colPercents",
                 "totPercents"),
               initialValues=c(
                 dialog.values$initial.rowPct,
                 dialog.values$initial.colPct,
                 dialog.values$initial.totPct),  ##c("0","0","0"),
               labels=gettextRcmdr(c(
                 "Row percentages",
                 "Column percentages",
                 "Total percentages")),
               title=gettextRcmdr("Compute Percentages"))


    checkBoxes(frame="testsFrame",
               boxes=c(
                 "chisqTest",
                 "chisqComponents",
                 "chiComp",
                 "expFreq",
                 "fisherTest"),
               initialValues=c(
                 dialog.values$initial.chisq,
                 dialog.values$initial.chisqComp,
                 dialog.values$initial.chiComp,
                 dialog.values$initial.expected,
                 dialog.values$initial.fisher), ## c("1", "1", "0", "0", "0"),
               labels=gettextRcmdr(c(
                 "Chi-square test of independence",
                 "Print chi-square components",
                 "Print chi components (residuals)",
                 "Print expected frequencies",
                 "Fisher's exact test")))


    tkgrid(getFrame(rowBox),
           labelRcmdr(variablesFrame, text="    "),
           getFrame(columnBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(top,
                   text=gettextRcmdr("Hypothesis Tests"),
                   fg="blue"),
           sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
  }

