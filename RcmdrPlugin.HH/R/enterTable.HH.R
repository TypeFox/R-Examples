"enterTable.HH" <-
function(){
    env <- environment()
    initializeDialog(title=gettextRcmdr("Enter Two-Way Table"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    setUpTable <- function(...){
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "tklabel(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep="")
            assign(col.varname, tclVar(j), envir=env)
            make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                    col.varname, ")", sep="")
            }
        eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows){   
            varname <- paste(".tab.", i, ".1", sep="") 
            assign(varname, tclVar("") , envir=env)
            row.varname <- paste(".rowname.", i, sep="")
            assign(row.varname, tclVar(i), envir=env)
            make.row <- paste("tkentry(.tableFrame, width='5', textvariable=",
                row.varname, ")", sep="")
            make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                varname, ")", sep="")
            for (j in 2:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(""), envir=env)
                make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='5', textvariable=", 
                    varname, ")", sep="")
                }
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
            }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
        }
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("2")
    rowsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=rowsValue,
        resolution=1, orient="horizontal", command=setUpTable)
    rowsShow <- tklabel(rowColFrame, textvariable=rowsValue, width=2, justify="right")
    colsValue <- tclVar("2")
    colsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=colsValue,
        resolution=1, orient="horizontal", command=setUpTable)
    colsShow <- tklabel(rowColFrame, textvariable=colsValue, width=2, justify="right")
    onOK <- function() {
        rowPct <- tclvalue(rowPercentsVariable)
        colPct <- tclvalue(colPercentsVariable)
        totPct <- tclvalue(totPercentsVariable)

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
            errorCondition(recall=enterTable.HH, message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
            return()
            }
        if (length(unique(row.names)) != nrows){
            errorCondition(recall=enterTable.HH, message=gettextRcmdr("Row names are not unique."))
            return()
            }     
        if (length(unique(col.names)) != ncols){
            errorCondition(recall=enterTable.HH, message=gettextRcmdr("Column names are not unique."))
            return()
            }     
##        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqCompVariable)
        chiComp <- tclvalue(chiCompVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
        closeDialog()

        command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
            ", byrow=TRUE)", sep="")
        justDoIt(paste(".Table <- ", command, sep=""))
        command <- paste("c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), ")", sep="")
        justDoIt(paste("rownames(.Table) <- ", command, sep=""))
        logger(paste("rownames(.Table) <- ", command, sep=""))
        command <- paste("c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), ")", sep="")
        justDoIt(paste("colnames(.Table) <- ", command, sep=""))
        logger(paste("colnames(.Table) <- ", command, sep=""))
        doItAndPrint(".Table  # Counts")

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
    
##     radioButtons(name="percents",
##                  buttons=c("rowPercents",
##                    "columnPercents",
##                    "nonePercents"),
##                  values=c("row", "column", "none"),
##                  initialValue="none",
##                  labels=gettextRcmdr(c
##                    ("Row percentages",
##                     "Column percentages",
##                     "No percentages")),
##                  title=gettextRcmdr("Compute Percentages"))
    
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

##     checkBoxes(frame="testsFrame",
##                boxes=c("chisq", "expFreq", "fisher"),
##                initialValues=c("1", "0", "0"),                               
##                labels=gettextRcmdr(c("Chi-square test of independence",
##                  "Print expected frequencies",
##                  "Fisher's exact test")))
    
    tkgrid(tklabel(rowColFrame,
                   text=gettextRcmdr("Number of Rows:")),
           rowsSlider, rowsShow, sticky="w")
    tkgrid(tklabel(rowColFrame,
                   text=gettextRcmdr("Number of Columns:")),
           colsSlider, colsShow, sticky="w")
    tkgrid(rowColFrame, sticky="w")
    tkgrid(tklabel(top,
                   text=gettextRcmdr("Enter counts:"),
                   fg="blue"), sticky="w")
    tkgrid(outerTableFrame, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Compute Percentages"), fg="blue"), sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(tklabel(top,
                   text=gettextRcmdr("Hypothesis Tests"),
                   fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=7, columns=2)
    }

