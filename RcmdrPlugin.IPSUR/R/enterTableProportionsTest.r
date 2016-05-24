# Last modified Feb 16, 2008


`enterTableMultiPropTest` <-
function () 
{
    require("abind")
    env <- environment()
    initializeDialog(title = gettextRcmdr("Enter table for multi-proportions test"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir = env)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir = env))
        assign(".tableFrame", tkframe(outerTableFrame), envir = env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "tklabel(.tableFrame, text='')"
        col.varname <- paste(".colname.", 1, sep = "")
        assign(col.varname, tclVar("Success"), envir = env)
        make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='10', textvariable=", 
            col.varname, ")", sep = "")
        col.varname <- paste(".colname.", 2, sep = "")
        assign(col.varname, tclVar("Failure"), envir = env)
        make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='10', textvariable=", 
            col.varname, ")", sep = "")
        eval(parse(text = paste("tkgrid(", make.col.names, ")", 
            sep = "")), envir = env)
        for (i in 1:nrows) {
            varname <- paste(".tab.", i, ".1", sep = "")
            assign(varname, tclVar(""), envir = env)
            row.varname <- paste(".rowname.", i, sep = "")
            assign(row.varname, tclVar(paste("Sample ", i, sep = "")), 
                envir = env)
            make.row <- paste("tkentry(.tableFrame, width='10', textvariable=", 
                row.varname, ")", sep = "")
            make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='10', textvariable=", 
                varname, ")", sep = "")
            for (j in 2:ncols) {
                varname <- paste(".tab.", i, ".", j, sep = "")
                assign(varname, tclVar(""), envir = env)
                make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='10', textvariable=", 
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
            errorCondition(recall = enterTableMultiPropTest, 
                message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), 
                  length(counts), nrows, ncols))
            return()
        }
        if (length(unique(row.names)) != nrows) {
            errorCondition(recall = enterTableMultiPropTest, 
                message = gettextRcmdr("Row names are not unique."))
            return()
        }
        if (length(unique(col.names)) != ncols) {
            errorCondition(recall = enterTableMultiPropTest, 
                message = gettextRcmdr("Column names are not unique."))
            return()
        }
        closeDialog()
        command <- paste(".Table <- matrix(c(", paste(counts, collapse = ","), 
            "), ", nrows, ", ", ncols, ", byrow=TRUE)", sep = "")
        justDoIt(command)
        logger(paste(command, sep = ""))
        command <- paste("c(", paste(paste("'", row.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("rownames(.Table) <- ", command, sep = ""))
        logger(paste("rownames(.Table) <- ", command, sep = ""))
        command <- paste("c(", paste(paste("'", col.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("colnames(.Table) <- ", command, sep = ""))
        logger(paste("colnames(.Table) <- ", command, sep = ""))
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceLevel)
        test <- as.character(tclvalue(testVariable))
        if (test == "normal") 
            doItAndPrint(paste("prop.test(.Table, alternative='", 
                alternative, "', conf.level=", level, ", correct=FALSE)", 
                sep = ""))
        else doItAndPrint(paste("prop.test(.Table, alternative='", 
            alternative, "', conf.level=", level, ", correct=TRUE)", 
            sep = ""))
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "prop.test")
    radioButtons(name = "alternative", buttons = c("twosided", 
        "less", "greater"), values = c("two.sided", "less", "greater"), 
        labels = gettextRcmdr(c("Two-sided", "Difference < 0 (samples=2)", 
            "Difference > 0 (samples=2)")), title = gettextRcmdr("Alternative Hypothesis"))
    confidenceFrame <- tkframe(top)
    confidenceLevel <- tclVar("0.95")
    confidenceField <- tkentry(confidenceFrame, width = "6", 
        textvariable = confidenceLevel)
    radioButtons(name = "test", buttons = c("normal", "corrected"), 
        labels = gettextRcmdr(c("Normal approximation", "Normal approximation with\ncontinuity correction (samples=2)")), 
        title = gettextRcmdr("Type of Test"))
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Rows (samples):")), 
        rowsSlider, rowsShow, sticky = "w")
    tkgrid(rowColFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter counts:"), 
        fg = "blue"), sticky = "w")
    tkgrid(outerTableFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("\nOptions:"), fg = "blue"), 
        sticky = "w")
    tkgrid(tklabel(confidenceFrame, text = gettextRcmdr("Confidence Level (samples=2): ")), 
        confidenceField, sticky = "w")
    tkgrid(confidenceFrame, sticky = "w")
    tkgrid(alternativeFrame, sticky = "nw")
    tkgrid(testFrame, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 5, columns = 2)
}



`enterTableSinglePropTest` <-
function () 
{
    require("abind")
    env <- environment()
    initializeDialog(title = gettextRcmdr("Enter table for single-proportion test"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir = env)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir = env))
        assign(".tableFrame", tkframe(outerTableFrame), envir = env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "tklabel(.tableFrame, text='')"
        col.varname <- paste(".colname.", 1, sep = "")
        assign(col.varname, tclVar("Success"), envir = env)
        make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='7', textvariable=", 
            col.varname, ")", sep = "")
        col.varname <- paste(".colname.", 2, sep = "")
        assign(col.varname, tclVar("Failure"), envir = env)
        make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='7', textvariable=", 
            col.varname, ")", sep = "")
        eval(parse(text = paste("tkgrid(", make.col.names, ")", 
            sep = "")), envir = env)
        for (i in 1:nrows) {
            varname <- paste(".tab.", i, ".1", sep = "")
            assign(varname, tclVar(""), envir = env)
            row.varname <- paste(".rowname.", i, sep = "")
            assign(row.varname, tclVar("Counts:"), envir = env)
            make.row <- paste("tkentry(.tableFrame, width='7', textvariable=", 
                row.varname, ")", sep = "")
            make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='7', textvariable=", 
                varname, ")", sep = "")
            for (j in 2:ncols) {
                varname <- paste(".tab.", i, ".", j, sep = "")
                assign(varname, tclVar(""), envir = env)
                make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='7', textvariable=", 
                  varname, ")", sep = "")
            }
            eval(parse(text = paste("tkgrid(", make.row, ")", 
                sep = "")), envir = env)
        }
        tkgrid(get(".tableFrame", envir = env), sticky = "w")
    }
    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("1")
    rowsSlider <- tkscale(rowColFrame, from = 1, to = 1, showvalue = FALSE, 
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
            errorCondition(recall = enterTableSinglePropTest, 
                message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), 
                  length(counts), nrows, ncols))
            return()
        }
        if (length(unique(col.names)) != ncols) {
            errorCondition(recall = enterTableSinglePropTest, 
                message = gettextRcmdr("Column names are not unique."))
            return()
        }
        closeDialog()
        command <- paste(".Table <- matrix(c(", paste(counts, collapse = ","), 
            "), ", nrows, ", ", ncols, ", byrow=TRUE)", sep = "")
        justDoIt(command)
        logger(paste(command, sep = ""))
        command <- paste("c(", paste(paste("'", row.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("rownames(.Table) <- ", command, sep = ""))
        logger(paste("rownames(.Table) <- ", command, sep = ""))
        command <- paste("c(", paste(paste("'", col.names, "'", 
            sep = ""), collapse = ", "), ")", sep = "")
        justDoIt(paste("colnames(.Table) <- ", command, sep = ""))
        logger(paste("colnames(.Table) <- ", command, sep = ""))
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceLevel)
        test <- as.character(tclvalue(testVariable))
        p <- tclvalue(pVariable)
        if (test == "normal") 
            doItAndPrint(paste("prop.test(rbind(.Table), alternative='", 
                alternative, "', p=", p, ", conf.level=", level, 
                ", correct=FALSE)", sep = ""))
        else if (test == "corrected") 
            doItAndPrint(paste("prop.test(rbind(.Table), alternative='", 
                alternative, "', p=", p, ", conf.level=", level, 
                ", correct=TRUE)", sep = ""))
        else doItAndPrint(paste("binom.test(rbind(.Table), alternative='", 
            alternative, "', p=", p, ", conf.level=", level, 
            ")", sep = ""))
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "prop.test")
    radioButtons(top, name = "alternative", buttons = c("twosided", 
        "less", "greater"), values = c("two.sided", "less", "greater"), 
        labels = gettextRcmdr(c("Population proportion = p0", 
            "Population proportion < p0", "Population proportion > p0")), 
        title = gettextRcmdr("Alternative Hypothesis"))
    rightFrame <- tkframe(top)
    confidenceFrame <- tkframe(top)
    confidenceLevel <- tclVar("0.95")
    confidenceField <- tkentry(confidenceFrame, width = "6", 
        textvariable = confidenceLevel)
    pFrame <- tkframe(top)
    pVariable <- tclVar("0.5")
    pField <- tkentry(pFrame, width = "6", textvariable = pVariable)
    radioButtons(name = "test", buttons = c("normal", "corrected", 
        "exact"), labels = gettextRcmdr(c("Normal approximation", 
        "Normal approximation with\ncontinuity correction", "Exact binomial")), 
        title = gettextRcmdr("Type of Test"))
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Rows (single sample):")), 
        rowsSlider, rowsShow, sticky = "w")
    tkgrid(rowColFrame, sticky = "w")
    tkgrid(outerTableFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("\nOptions:"), fg = "blue"), 
        sticky = "w")
    tkgrid(tklabel(pFrame, text = gettextRcmdr("Null hypothesis: p0 = ")), 
        pField, sticky = "w")
    tkgrid(pFrame, sticky = "w")
    tkgrid(tklabel(confidenceFrame, text = gettextRcmdr("Confidence Level:     ")), 
        confidenceField, sticky = "w")
    tkgrid(confidenceFrame, sticky = "w")
    tkgrid(alternativeFrame, sticky = "nw")
    tkgrid(testFrame, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 5, columns = 2)
}
