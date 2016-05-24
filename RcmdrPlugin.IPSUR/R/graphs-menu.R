# Last modified Feb 16, 2008


`barGraph.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Bar Graph"))
    variableBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable (pick one)"))
    .groups <- FALSE
    onOK <- function() {
        variable <- getSelection(variableBox)
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = barGraph.ipsur, message = gettextRcmdr("You must select a variable"))
            return()
        }
        title <- tclvalue(titleVariable)
        leg <- tclvalue(legendVariable) == "1"
        prop <- tclvalue(propVariable) == "1"
        besid <- tclvalue(typeVariable) == "beside"
        if (tclvalue(coloVariable) == 1) {
            if (.groups == FALSE) {
                colcomm <- paste(", col=rainbow(length(table(", 
                  ActiveDataSet(), "$", variable, "))))", sep = "")
            }
            else {
                colcomm <- paste(", col=rainbow(length(table(", 
                  ActiveDataSet(), "$", .groups, "))))", sep = "")
            }
        }
        else {
            colcomm <- ", col=NULL)"
        }
        if (prop) {
            if (.groups == FALSE) {
                command <- paste("barplot(prop.table(table(", 
                  ActiveDataSet(), "$", variable, ")), main=\"", 
                  title, "\", xlab=\"", variable, "\", ylab=\"Relative Frequency\", legend.text=", 
                  leg, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
            else {
                command <- paste("barplot(prop.table(table(", 
                  ActiveDataSet(), "$", .groups, ", ", ActiveDataSet(), 
                  "$", variable, ")), main=\"", title, "\", xlab=\"", 
                  variable, "\", ylab=\"Relative Frequency\", legend.text=", 
                  leg, ", beside=", besid, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
        }
        else {
            if (.groups == FALSE) {
                command <- paste("barplot(table(", ActiveDataSet(), 
                  "$", variable, "), main=\"", title, "\", xlab=\"", 
                  variable, "\", ylab=\"Frequency\", legend.text=", 
                  leg, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
            else {
                command <- paste("barplot(table(", ActiveDataSet(), 
                  "$", .groups, ", ", ActiveDataSet(), "$", variable, 
                  "), main=\"", title, "\", xlab=\"", variable, 
                  "\", ylab=\"Frequency\", legend.text=", leg, 
                  ", beside=", besid, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(barGraph.ipsur)
    OKCancelHelp(helpSubject = "barplot")
    radioButtons(name = "type", buttons = c("segmented", "beside"), 
        labels = gettextRcmdr(c("Stacked bars", "Side-by-side bars")), 
        title = gettextRcmdr("Display groups with:"))
    optionsFrame <- tkframe(top)
    legendVariable <- tclVar("1")
    legendCheckBox <- tkcheckbutton(optionsFrame, variable = legendVariable)
    propVariable <- tclVar("0")
    propCheckBox <- tkcheckbutton(optionsFrame, variable = propVariable)
    coloVariable <- tclVar("0")
    coloCheckBox <- tkcheckbutton(optionsFrame, variable = coloVariable)
    titleFrame <- tkframe(top)
    titleVariable <- tclVar(gettextRcmdr(""))
    titleField <- tkentry(titleFrame, width = "40", textvariable = titleVariable)
    tkgrid(tklabel(titleFrame, text = gettextRcmdr("Title: "), 
        fg = "blue"), titleField, sticky = "w")
    tkgrid(titleFrame, sticky = "w")
    tkgrid(getFrame(variableBox), sticky = "nw")
    tkgrid(tklabel(top, text = gettextRcmdr("Options:"), fg = "blue"), 
        sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Relative Frequencies: "), 
        justify = "left"), propCheckBox, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Rainbow: "), 
        justify = "left"), coloCheckBox, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Legend: "), 
        justify = "left"), legendCheckBox, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(typeFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}




`barPlotSumTable` <-
function () 
{
    require("abind")
    env <- environment()
    initializeDialog(title = gettextRcmdr("Bar Graph for Summarized Data"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir = env)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir = env))
        assign(".tableFrame", tkframe(outerTableFrame), envir = env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- as.numeric(tclvalue(colsValue))
        make.col.names <- "tklabel(.tableFrame, text='')"
        for (j in 1:ncols) {
            col.varname <- paste(".colname.", j, sep = "")
            assign(col.varname, tclVar(paste("Response ", j, 
                sep = "")), envir = env)
            make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='10', textvariable=", 
                col.varname, ")", sep = "")
        }
        eval(parse(text = paste("tkgrid(", make.col.names, ")", 
            sep = "")), envir = env)
        for (i in 1:nrows) {
            varname <- paste(".tab.", i, ".1", sep = "")
            assign(varname, tclVar(""), envir = env)
            row.varname <- paste(".rowname.", i, sep = "")
            assign(row.varname, tclVar(paste("Group ", i, sep = "")), 
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
    rowsSlider <- tkscale(rowColFrame, from = 1, to = 10, showvalue = FALSE, 
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
            errorCondition(recall = barPlotSumTable, message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), 
                length(counts), nrows, ncols))
            return()
        }
        if (length(unique(row.names)) != nrows) {
            errorCondition(recall = barPlotSumTable, message = gettextRcmdr("Row names are not unique."))
            return()
        }
        if (length(unique(col.names)) != ncols) {
            errorCondition(recall = barPlotSumTable, message = gettextRcmdr("Column names are not unique."))
            return()
        }
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
        title <- tclvalue(titleVariable)
        leg <- tclvalue(legendVariable) == "1"
        prop <- tclvalue(propVariable) == "1"
        besid <- tclvalue(typeVariable) == "beside"
        if (tclvalue(coloVariable) == 1) {
            colcomm <- paste(", col=rainbow(", nrows, "))", sep = "")
        }
        else {
            colcomm <- ", col=NULL)"
        }
        if (prop) {
            if (nrows == 1) {
                command <- paste("barplot(prop.table(.Table), main=\"", 
                  title, "\", ylab=\"Relative Frequency\"", colcomm, 
                  sep = "")
                logger(command)
                justDoIt(command)
            }
            else {
                command <- paste("barplot(prop.table(.Table), main=\"", 
                  title, "\", ylab=\"Relative Frequency\", legend.text=", 
                  leg, ", beside=", besid, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
        }
        else {
            if (nrows == 1) {
                command <- paste("barplot(.Table, main=\"", title, 
                  "\", ylab=\"Frequency\"", colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
            else {
                command <- paste("barplot(.Table, main=\"", title, 
                  "\", ylab=\"Frequency\", legend.text=", leg, 
                  ", beside=", besid, colcomm, sep = "")
                logger(command)
                justDoIt(command)
            }
        }
    }
    OKCancelHelp(helpSubject = "barplot")
    titleFrame <- tkframe(top)
    titleVariable <- tclVar(gettextRcmdr(""))
    titleField <- tkentry(titleFrame, width = "40", textvariable = titleVariable)
    tkgrid(tklabel(titleFrame, text = gettextRcmdr("Title: "), 
        fg = "blue"), titleField, sticky = "w")
    tkgrid(titleFrame, sticky = "w")
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Columns (reponses):")), 
        colsSlider, colsShow, sticky = "w")
    tkgrid(tklabel(rowColFrame, text = gettextRcmdr("Number of Rows (groups):")), 
        rowsSlider, rowsShow, sticky = "w")
    tkgrid(rowColFrame, sticky = "w")
    tkgrid(tklabel(top, text = gettextRcmdr("Enter counts:"), 
        fg = "blue"), sticky = "w")
    tkgrid(outerTableFrame, sticky = "w")
    radioButtons(name = "type", buttons = c("segmented", "beside"), 
        labels = gettextRcmdr(c("Stacked bars", "Side-by-side bars")), 
        title = gettextRcmdr("Display groups with:"))
    optionsFrame <- tkframe(top)
    legendVariable <- tclVar("1")
    legendCheckBox <- tkcheckbutton(optionsFrame, variable = legendVariable)
    propVariable <- tclVar("0")
    propCheckBox <- tkcheckbutton(optionsFrame, variable = propVariable)
    coloVariable <- tclVar("0")
    coloCheckBox <- tkcheckbutton(optionsFrame, variable = coloVariable)
    tkgrid(tklabel(top, text = gettextRcmdr("Options:"), fg = "blue"), 
        sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Relative Frequencies: "), 
        justify = "left"), propCheckBox, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Rainbow: "), 
        justify = "left"), coloCheckBox, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Legend (groups>1): "), 
        justify = "left"), legendCheckBox, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(typeFrame, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 7, columns = 2)
}




`boxPlot.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Boxplot"))
    xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"))
    identifyVariable <- tclVar("0")
    identifyFrame <- tkframe(top)
    identifyCheckBox <- tkcheckbutton(identifyFrame, variable = identifyVariable)
    .groups <- FALSE
    onOK <- function() {
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = boxPlot.ipsur, message = gettextRcmdr("You must select a variable"))
            return()
        }
        identifyPoints <- "1" == tclvalue(identifyVariable)
        .activeDataSet <- ActiveDataSet()
        var <- paste(.activeDataSet, "$", x, sep = "")
        title <- tclvalue(titleVariable)
        horz <- tclvalue(horzOptionVariable)
        notch <- tclvalue(notchOptionVariable)
        varwid <- tclvalue(varwidthOptionVariable)
        hflag <- (horz == 1)
        nflag <- (notch == 1)
        vwflag <- (varwid == 1)
        if (.groups == FALSE) {
            if (hflag) {
                command <- (paste("boxplot(", var, ", xlab=\"", 
                  x, "\", main=\"", title, "\", notch=", nflag, 
                  ", warwidth=", vwflag, ", \n                              horizontal=", 
                  hflag, ")", sep = ""))
            }
            else {
                command <- (paste("boxplot(", var, ", ylab=\"", 
                  x, "\", main=\"", title, "\", notch=", nflag, 
                  ", warwidth=", vwflag, ", \n                              horizontal=", 
                  hflag, ")", sep = ""))
            }
            logger(command)
            justDoIt(command)
            if (identifyPoints) {
                RcmdrTkmessageBox(title = "Identify Points", 
                  message = gettextRcmdr("Use left mouse button to identify points,\nright button to exit."), 
                  icon = "info", type = "ok")
                doItAndPrint(paste("identify(rep(1, length(", 
                  var, ")), ", var, ", rownames(", .activeDataSet, 
                  "))", sep = ""))
            }
        }
        else {
            if (hflag) {
                command <- (paste("boxplot(", x, "~", .groups, 
                  ", xlab=\"", x, "\", xlab=\"", .groups, "\"", 
                  ", main=\"", title, "\", notch=", nflag, ", varwidth=", 
                  vwflag, ", horizontal=", hflag, ", data=", 
                  .activeDataSet, ")", sep = ""))
            }
            else {
                command <- (paste("boxplot(", x, "~", .groups, 
                  ", ylab=\"", x, "\", xlab=\"", .groups, "\"", 
                  ", main=\"", title, "\", notch=", nflag, ", varwidth=", 
                  vwflag, ", horizontal=", hflag, ", data=", 
                  .activeDataSet, ")", sep = ""))
            }
            logger(command)
            justDoIt(command)
            if (identifyPoints) {
                RcmdrTkmessageBox(title = "Identify Points", 
                  message = gettextRcmdr("Use left mouse button to identify points,\nright button to exit."), 
                  icon = "info", type = "ok")
                doItAndPrint(paste("identify(", .activeDataSet, 
                  "$", .groups, ", ", var, ", rownames(", .activeDataSet, 
                  "))", sep = ""))
            }
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(boxPlot.ipsur)
    OKCancelHelp(helpSubject = "boxplot")
    optionsFrame <- tkframe(top)
    checkBoxes(frame = "optionsFrame", boxes = c("horzOption", 
        "notchOption", "varwidthOption"), initialValues = c("1", 
        "0", "1"), labels = gettextRcmdr(c("Horizontal", "Notches", 
        "Variable Box Width")))
    titleFrame <- tkframe(top)
    titleVariable <- tclVar(gettextRcmdr(""))
    titleField <- tkentry(titleFrame, width = "40", textvariable = titleVariable)
    tkgrid(tklabel(titleFrame, text = gettextRcmdr("Title: "), 
        fg = "blue"), titleField, sticky = "w")
    tkgrid(titleFrame, sticky = "w")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(tklabel(top, text = gettextRcmdr("Options"), fg = "blue"), 
        sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(tklabel(identifyFrame, text = gettextRcmdr("Identify outliers with mouse"), 
        justify = "left"), identifyCheckBox, sticky = "w")
    tkgrid(identifyFrame, stick = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 4, columns = 1)
}




`paretoChart` <-
function () 
{
    require("qcc")
    initializeDialog(title = gettextRcmdr("Pareto Chart"))
    variableBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable (pick one)"))
    onOK <- function() {
        variable <- getSelection(variableBox)
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = paretoChart, message = gettextRcmdr("You must select a variable"))
            return()
        }
        title <- tclvalue(titleVariable)
        colo <- tclvalue(coloOptionVariable)
        if (colo == 1) {
            colcomm <- paste(", col=rainbow(length(table(", ActiveDataSet(), 
                "$", variable, "))))", sep = "")
        }
        else {
            colcomm <- ", col=8)"
        }
        command <- paste("pareto.chart(table(", ActiveDataSet(), 
            "$", variable, "), main=\"", title, "\", ylab=\"Frequency\"", 
            colcomm, sep = "")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pareto.chart")
    optionsFrame <- tkframe(top)
    checkBoxes(frame = "optionsFrame", boxes = c("coloOption"), 
        initialValues = c("0"), labels = gettextRcmdr(c("Rainbow")))
    titleFrame <- tkframe(top)
    titleVariable <- tclVar(gettextRcmdr(""))
    titleField <- tkentry(titleFrame, width = "40", textvariable = titleVariable)
    tkgrid(tklabel(titleFrame, text = gettextRcmdr("Title: "), 
        fg = "blue"), titleField, sticky = "w")
    tkgrid(titleFrame, sticky = "w")
    tkgrid(getFrame(variableBox), sticky = "nw")
    tkgrid(tklabel(top, text = gettextRcmdr("Options"), fg = "blue"), 
        sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}



`stripChart` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Strip chart"))
    .numeric <- Numeric()
    optionsFrame <- tkframe(top)
    groupsFrame <- tkframe(top)
    xBox <- variableListBox(top, .numeric, title = gettextRcmdr("Variable (pick one)"))
    slider1Value <- tclVar("1")
    slider1 <- tkscale(optionsFrame, from = 1, to = 25, showvalue = TRUE, 
        variable = slider1Value, resolution = 1, orient = "horizontal")
    slider2Value <- tclVar("1")
    slider2 <- tkscale(optionsFrame, from = 0.1, to = 2.5, showvalue = TRUE, 
        variable = slider2Value, resolution = 0.1, orient = "horizontal")
    .groups <- FALSE
    onOK <- function() {
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = stripChart, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        title <- tclvalue(titleVariable)
        method <- tclvalue(methodVariable)
        jitter <- if (method == "jitter") 
            paste(", jitter=", tclvalue(jitterVariable), sep = "")
        else ""
        offset <- if (method == "stack") 
            paste(", offset=", tclvalue(offsetVariable), sep = "")
        else ""
        pch <- as.numeric(tclvalue(slider1Value))
        cex <- as.numeric(tclvalue(slider2Value))
        tkdestroy(top)
        if (.groups == FALSE) {
            doItAndPrint(paste("stripchart(", .activeDataSet, 
                "$", x, ", method=\"", method, "\"", jitter, 
                offset, ", main=\"", title, "\", pch=", pch, 
                ", cex=", cex, ")", sep = ""))
        }
        else {
            doItAndPrint(paste("stripchart(", .activeDataSet, 
                "$", x, "~", .activeDataSet, "$", .groups, ", method=\"", 
                method, "\"", jitter, offset, ", main=\"", title, 
                "\", pch=", pch, ", cex=", cex, ")", sep = ""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(stripChart)
    OKCancelHelp(helpSubject = "stripchart")
    methodVariable <- tclVar("overplot")
    overplotButton <- tkradiobutton(optionsFrame, variable = methodVariable, 
        value = "overplot")
    jitterButton <- tkradiobutton(optionsFrame, variable = methodVariable, 
        value = "jitter")
    stackButton <- tkradiobutton(optionsFrame, variable = methodVariable, 
        value = "stack")
    jitterVariable <- tclVar("0.1")
    offsetVariable <- tclVar("1/3")
    jitterEntry <- tkentry(optionsFrame, width = "6", textvariable = jitterVariable)
    offsetEntry <- tkentry(optionsFrame, width = "6", textvariable = offsetVariable)
    titleFrame <- tkframe(top)
    titleVariable <- tclVar(gettextRcmdr(""))
    titleField <- tkentry(titleFrame, width = "30", textvariable = titleVariable)
    tkgrid(tklabel(titleFrame, text = gettextRcmdr("Title: "), 
        fg = "blue"), titleField, sticky = "w")
    tkgrid(titleFrame, sticky = "w")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Plot method:"), 
        fg = "blue"), sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Overplot")), 
        overplotButton, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Jitter")), 
        jitterButton, tklabel(optionsFrame, text = gettextRcmdr("   amount:")), 
        jitterEntry, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Stack")), 
        stackButton, tklabel(optionsFrame, text = gettextRcmdr("   offset:")), 
        offsetEntry, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Plot Character:"), 
        fg = "blue"), tklabel(optionsFrame, text = gettextRcmdr("")), 
        slider1, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Character Size:"), 
        fg = "blue"), tklabel(optionsFrame, text = gettextRcmdr("")), 
        slider2, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}
