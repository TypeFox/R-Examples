# Last modified Feb 16, 2008

`enterTable.ipsur` <-
function () 
{
    require("abind")
    env <- environment()
    initializeDialog(title = gettextRcmdr("Enter Two-Way Table"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir = env)
    simulateFrame <- tkframe(top)
    simulateVariable <- tclVar("0")
    simulateCheckBox <- tkcheckbutton(simulateFrame, variable = simulateVariable)
    simulate <- tclVar("2000")
    simulateEntry <- tkentry(simulateFrame, width = "6", textvariable = simulate)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir = env))
        assign(".tableFrame", tkframe(outerTableFrame), envir = env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "tklabel(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep = "")
            assign(col.varname, tclVar(j), envir = env)
            make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                col.varname, ")", sep = "")
        }
        eval(parse(text = paste("tkgrid(", make.col.names, ")", 
            sep = "")), envir = env)
        for (i in 1:nrows) {
            varname <- paste(".tab.", i, ".1", sep = "")
            assign(varname, tclVar(""), envir = env)
            row.varname <- paste(".rowname.", i, sep = "")
            assign(row.varname, tclVar(i), envir = env)
            make.row <- paste("tkentry(.tableFrame, width='5', textvariable=", 
                row.varname, ")", sep = "")
            make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                varname, ")", sep = "")
            for (j in 2:ncols) {
                varname <- paste(".tab.", i, ".", j, sep = "")
                assign(varname, tclVar(""), envir = env)
                make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                  varname, ")", sep = "")
            }
            eval(parse(text = paste("tkgrid(", make.row, ")", 
                sep = "")), envir = env)
        }
        tkgrid(get(".tableFrame", envir = env), sticky = "w")
    }
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("2")
    rowsSlider <- tkscale(rowColFrame, from = 2, to = 10, showvalue = FALSE, 
        variable = rowsValue, resolution = 1, orient = "horizontal", 
        command = setUpTable)
    rowsShow <- tklabel(rowColFrame, textvariable = rowsValue, 
        width = 2, justify = "right")
    colsValue <- tclVar("2")
    colsSlider <- tkscale(rowColFrame, from = 2, to = 10, showvalue = FALSE, 
        variable = colsValue, resolution = 1, orient = "horizontal", 
        command = setUpTable)
    colsShow <- tklabel(rowColFrame, textvariable = colsValue, 
        width = 2, justify = "right")
    onOK <- function() {
        margins <- tclvalue(marginsVariable)
        sims <- tclvalue(simulateVariable)
        B <- tclvalue(simulate)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        cell <- 0
        counts <- rep(NA, nrows * ncols)
        row.names <- rep("", nrows)
        col.names <- rep("", ncols)
        for (i in 1:nrows) row.names[i] <- eval(parse(text = paste("tclvalue(", 
            paste(".rowname.", i, sep = ""), ")", sep = "")))
        for (j in 1:ncols) col.names[j] <- eval(parse(text = paste("tclvalue(", 
            paste(".colname.", j, sep = ""), ")", sep = "")))
        for (i in 1:nrows) {
            for (j in 1:ncols) {
                cell <- cell + 1
                varname <- paste(".tab.", i, ".", j, sep = "")
                counts[cell] <- as.numeric(eval(parse(text = paste("tclvalue(", 
                  varname, ")", sep = ""))))
            }
        }
        counts <- na.omit(counts)
        if (length(counts) != nrows * ncols) {
            errorCondition(recall = enterTable.ipsur, message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), 
                length(counts), nrows, ncols))
            return()
        }
        if (length(unique(row.names)) != nrows) {
            errorCondition(recall = enterTable.ipsur, message = gettextRcmdr("Row names are not unique."))
            return()
        }
        if (length(unique(col.names)) != ncols) {
            errorCondition(recall = enterTable.ipsur, message = gettextRcmdr("Column names are not unique."))
            return()
        }
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherVariable)
        closeDialog()
        command <- paste(".Table <- matrix(c(", paste(counts, collapse = ","), 
            "), ", nrows, ", ", ncols, ", byrow=TRUE)", sep = "")
        justDoIt(command)
        logger(paste(".Table <- ", command, sep = ""))
        command <- paste("c(", paste(paste("'", row.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("rownames(.Table) <- ", command, sep = ""))
        logger(paste("rownames(.Table) <- ", command, sep = ""))
        command <- paste("c(", paste(paste("'", col.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("colnames(.Table) <- ", command, sep = ""))
        logger(paste("colnames(.Table) <- ", command, sep = ""))
        if (margins == 1) {
            doItAndPrint("addmargins(.Table) # Counts with Marginal Distributions")
        }
        else {
            doItAndPrint(".Table  # Counts")
        }
        if (percents == "row") 
            doItAndPrint("rowPercents(.Table) # Row Percentages")
        if (percents == "column") 
            doItAndPrint("colPercents(.Table) # Column Percentages")
        if (percents == "total") 
            doItAndPrint("totPercents(.Table) # Percentage of Total")
        if (chisq == 1) {
            if (sims == 0) {
                command <- ".Test <- chisq.test(.Table, correct=FALSE)"
            }
            else {
                command <- paste(".Test <- chisq.test(.Table, correct=FALSE, simulate.p.value=TRUE, B=", 
                  B, ")", sep = "")
            }
            logger(paste(".Test <- ", command, sep = ""))
            justDoIt(command)
            doItAndPrint(".Test")
            if (expected == 1) 
                doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) 
                warnText <- paste(nlt1, gettextRcmdr("expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) 
                warnText <- paste(warnText, "\n", nlt5, gettextRcmdr(" expected frequencies are less than 5"), 
                  sep = "")
            if (!is.null(warnText)) 
                Message(message = warnText, type = "warning")
            if (chisqComp == 1) {
                command <- "round(.Test$residuals^2, 2) # Chi-square Components"
                doItAndPrint(command)
            }
            logger("remove(.Test)")
            remove(.Test, envir = .GlobalEnv)
        }
        if (fisher == 1) 
            doItAndPrint("fisher.test(.Table)")
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "chisq.test")
    radioButtons(name = "percents", buttons = c("rowPercents", 
        "columnPercents", "totalPercents", "nonePercents"), values = c("row", 
        "column", "total", "none"), initialValue = "none", labels = gettextRcmdr(c("Row percentages", 
        "Column percentages", "Percentages of total", "No percentages")), 
        title = gettextRcmdr("Compute Percentages"))
    checkBoxes(frame = "testsFrame", boxes = c("chisq", "chisqComponents", 
        "expFreq", "fisher"), initialValues = c("1", "0", "0", 
        "0"), labels = gettextRcmdr(c("Chi-square test of independence", 
        "Components of chi-square statistic", "Print expected frequencies", 
        "Fisher's exact test")))
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Rows:")), 
        rowsSlider, rowsShow, sticky = "w")
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Columns:")), 
        colsSlider, colsShow, sticky = "w")
    tkgrid(rowColFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter counts:"), 
        fg = "blue"), sticky = "w")
    tkgrid(outerTableFrame, sticky = "w")
    checkBoxes(frame = "marginsFrame", boxes = c("margins"), 
        initialValues = c("1"), labels = gettextRcmdr(c("Add Marginal Distributions")))
    tkgrid(marginsFrame, sticky = "w")
    tkgrid(percentsFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Hypothesis Tests"), 
        fg = "blue"), sticky = "w")
    tkgrid(testsFrame, sticky = "w")
    tkgrid(tklabel(simulateFrame, text = gettextRcmdr("Simulate p-value")), 
        simulateCheckBox, tklabel(simulateFrame, text = gettextRcmdr(" Iterations:")), 
        simulateEntry, sticky = "w")
    tkgrid(simulateFrame, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 7, columns = 2)
}


`twoWayTable.ipsur` <-
function () 
{
    require("abind")
    initializeDialog(title = gettextRcmdr("Two-Way Table"))
    variablesFrame <- tkframe(top)
    simulateFrame <- tkframe(top)
    simulateVariable <- tclVar("0")
    simulateCheckBox <- tkcheckbutton(simulateFrame, variable = simulateVariable)
    simulate <- tclVar("2000")
    simulateEntry <- tkentry(simulateFrame, width = "6", textvariable = simulate)
    .factors <- Factors()
    rowBox <- variableListBox(variablesFrame, .factors, title = gettextRcmdr("Row variable (pick one)"))
    columnBox <- variableListBox(variablesFrame, .factors, title = gettextRcmdr("Column variable (pick one)"))
    subsetBox()
    onOK <- function() {
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
        if (length(row) == 0 || length(column) == 0) {
            errorCondition(recall = twoWayTable.ipsur, message = gettextRcmdr("You must select two variables."))
            return()
        }
        if (row == column) {
            errorCondition(recall = twoWayTable.ipsur, message = gettextRcmdr("Row and column variables are the same."))
            return()
        }
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
        margins <- tclvalue(marginsVariable)
        sims <- tclvalue(simulateVariable)
        B <- tclvalue(simulate)
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) 
            ""
        else paste(", subset=", subset, sep = "")
        closeDialog()
        command <- paste(".Table <- xtabs(~", row, "+", column, ", data=", 
            ActiveDataSet(), subset, ")", sep = "")
        logger(paste(".Table <- ", command, sep = ""))
        justDoIt(command)
        if (margins == 1) {
            doItAndPrint("addmargins(.Table) # Table with Marginal Distributions")
        }
        else {
            doItAndPrint(".Table")
        }
        if (percents == "row") 
            doItAndPrint("rowPercents(.Table) # Row Percentages")
        if (percents == "column") 
            doItAndPrint("colPercents(.Table) # Column Percentages")
        if (percents == "total") 
            doItAndPrint("totPercents(.Table) # Percentage of Total")
        if (chisq == 1) {
            if (sims == 0) {
                command <- ".Test <- chisq.test(.Table, correct=FALSE)"
            }
            else {
                command <- paste(".Test <- chisq.test(.Table, correct=FALSE, simulate.p.value=TRUE, B=", 
                  B, ")", sep = "")
            }
            logger(paste(".Test <- ", command, sep = ""))
            justDoIt(command)
            doItAndPrint(".Test")
            if (expected == 1) 
                doItAndPrint(".Test$expected # Expected Counts")
            warnText <- NULL
            if (0 < (nlt1 <- sum(.Test$expected < 1))) 
                warnText <- paste(nlt1, gettextRcmdr("expected frequencies are less than 1"))
            if (0 < (nlt5 <- sum(.Test$expected < 5))) 
                warnText <- paste(warnText, "\n", nlt5, gettextRcmdr(" expected frequencies are less than 5"), 
                  sep = "")
            if (!is.null(warnText)) 
                Message(message = warnText, type = "warning")
            if (chisqComp == 1) {
                command <- "round(.Test$residuals^2, 2) # Chi-square Components"
                doItAndPrint(command)
            }
            logger("remove(.Test)")
            remove(.Test, envir = .GlobalEnv)
        }
        if (fisher == 1) 
            doItAndPrint("fisher.test(.Table)")
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "xtabs")
    radioButtons(name = "percents", buttons = c("rowPercents", 
        "columnPercents", "totalPercents", "nonePercents"), values = c("row", 
        "column", "total", "none"), initialValue = "none", labels = gettextRcmdr(c("Row percentages", 
        "Column percentages", "Percentages of total", "No percentages")), 
        title = gettextRcmdr("Compute Percentages"))
    checkBoxes(frame = "testsFrame", boxes = c("chisqTest", "chisqComponents", 
        "expFreq", "fisherTest"), initialValues = c("1", "0", 
        "0", "0"), labels = gettextRcmdr(c("Chi-square test of independence", 
        "Components of chi-square statistic", "Print expected frequencies", 
        "Fisher's exact test")))
    tkgrid(getFrame(rowBox), tklabel(variablesFrame, text = "    "), 
        getFrame(columnBox), sticky = "nw")
    tkgrid(variablesFrame, sticky = "w")
    checkBoxes(frame = "marginsFrame", boxes = c("margins"), 
        initialValues = c("1"), labels = gettextRcmdr(c("Add Marginal Distributions")))
    tkgrid(marginsFrame, sticky = "w")
    tkgrid(percentsFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Hypothesis Tests"), 
        fg = "blue"), sticky = "w")
    tkgrid(testsFrame, sticky = "w")
    tkgrid(tklabel(simulateFrame, text = gettextRcmdr("Simulate p-value")), 
        simulateCheckBox, tklabel(simulateFrame, text = gettextRcmdr(" Iterations:")), 
        simulateEntry, sticky = "w")
    tkgrid(simulateFrame, sticky = "w")
    tkgrid(subsetFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 6, columns = 1)
}
