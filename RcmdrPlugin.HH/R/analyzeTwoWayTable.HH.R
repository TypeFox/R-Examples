"analyzeTwoWayTable.HH" <-
function(){
    initializeDialog(title=gettextRcmdr("Analyze Two-Way Table"))
    variablesFrame <- tkframe(top)
    onOK <- function() {
        rowPct <- tclvalue(rowPercentsVariable)
        colPct <- tclvalue(colPercentsVariable)
        totPct <- tclvalue(totPercentsVariable)

        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqCompVariable)
        chiComp <- tclvalue(chiCompVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
        closeDialog()

        ## .Table <- data.matrix(ActiveDataSet())
        command <- ActiveDataSet()
        justDoIt(paste(".Table <- ", command, sep=""))
        doItAndPrint(".Table")

        if (rowPct == 1) doItAndPrint("rowPercents(.Table) # Row Percentages")
        if (colPct == 1) doItAndPrint("colPercents(.Table) # Column Percentages")        
        if (totPct == 1) doItAndPrint("totPercents(.Table) # Total Percentages")        

        if (chisq == 1) {
            command <- "chisq.test(.Table, correct=FALSE)"
            justDoIt(paste(".Test <- ", command, sep=""))
            doItAndPrint(".Test")
            if (chisqComp == 1) doItAndPrint("round(.Test$residuals^2, 2) # Chi-square Components")
            if (chiComp == 1) doItAndPrint("round(.Test$residuals, 2) # Chi Components (residuals)")
            if (expected == 1) doItAndPrint("round(.Test$expected, 2) # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
                gettextRcmdr("expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                gettextRcmdr(" expected frequencies are less than 5"), sep="")
            if (!is.null(warnText)) Message(message=warnText,
                type="warning")
            logger("remove(.Test)") 
            remove(.Test, envir=.GlobalEnv) 
            }
        if (fisher == 1) doItAndPrint("fisher.test(.Table)")

        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="chisq.test")

    checkBoxes(frame="percentsFrame",
               boxes=c("rowPercents",
                 "colPercents",
                 "totPercents"), 
               initialValues=c("0","0","0"), 
               labels=gettextRcmdr(c(
                 "Row percentages",
                 "Column percentages",
                 "Total percentages")))

    checkBoxes(frame="testsFrame",
               boxes=c(
                 "chisqTest",
                 "chisqComp",
                 "chiComp",
                 "expFreq",
                 "fisherTest"),
               initialValues=c("1", "1", "0", "0", "0"),
               labels=gettextRcmdr(c(
                 "Chi-square test of independence",
                 "Print chi-square components",
                 "Print chi components (residuals)",
                 "Print expected frequencies",
                 "Fisher's exact test")))

    tkgrid(variablesFrame, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Compute Percentages"), fg="blue"), sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
  }

