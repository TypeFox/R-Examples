varTableDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Variable One-Way Table"))

    vars <- colnames(meta(corpus))
    varBox <- variableListBox(top, vars,
                               title=.gettext("Variable:"),
                               initialSelection=0)

    radioButtons(name="what",
                 buttons=c("percent", "absolute"),
                 labels=c(.gettext("Percent"),
                          .gettext("Absolute counts")),
                 title=.gettext("Measure:"),
                 right.buttons=FALSE)

    plotFrame <- tkframe(top)

    tclPlotVar <- tclVar(1)
    plotButton <- tkcheckbutton(plotFrame, text=.gettext("Draw plot"), variable=tclPlotVar)

    tclVertVar <- tclVar(0)
    vertButton <- tkcheckbutton(plotFrame, text=.gettext("Vertical bars"), variable=tclVertVar)

    tclTitle <- tclVar(.gettext("Distribution of documents by %V"))
    titleEntry <- ttkentry(plotFrame, width=40, textvariable=tclTitle)

    onOK <- function() {
        var <- getSelection(varBox)
        what <- tclvalue(whatVariable)
        plot <- tclvalue(tclPlotVar) == 1
        vert <- tclvalue(tclVertVar) == 1
        title <- tclvalue(tclTitle)

        closeDialog()

        title <- gsub("%V", tolower(var), title)

        doItAndPrint(sprintf('absVarFreqs <- table(meta(corpus, "%s"), dnn="%s")', var, var))

        if(what == "percent") {
            doItAndPrint("varFreqs <- prop.table(absVarFreqs) * 100")
            ylab <- .gettext("% of documents")
        }
        else {
            doItAndPrint("varFreqs <- absVarFreqs")
            ylab <- .gettext("Number of documents")
        }

        if(plot) {
            if(vert)
                doItAndPrint(sprintf('barchart(varFreqs, horizontal=FALSE, scales=list(rot=90), ylab="%s", main="%s", auto.key=TRUE)',
                                     ylab, title))
            else
                doItAndPrint(sprintf('barchart(varFreqs, xlab="%s", main="%s", auto.key=TRUE)',
                                     ylab, title))
        }

        doItAndPrint("varFreqs <- addmargins(varFreqs)")

        if(what == "percent")
            setLastTable("varFreqs", paste(title, "(%)"))
        else
            setLastTable("varFreqs", title)

        doItAndPrint("varFreqs")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="varTableDlg")
    tkgrid(getFrame(varBox), sticky="w", pady=6, columnspan=2)
    tkgrid(whatFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(.titleLabel(plotFrame, text=.gettext("Plot:")),
           sticky="w", columnspan=2)
    tkgrid(plotButton, sticky="w", columnspan=2)
    tkgrid(vertButton, sticky="w", columnspan=2)
    tkgrid(labelRcmdr(plotFrame, text=.gettext("Title:")), titleEntry, sticky="w")
    tkgrid(plotFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=varBox$listbox)
}

varCrossTableDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }
    else if(nVars == 1) {
        .Message(message=.gettext("Corpus has only one variable."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Variable Two-Way Table"))

    vars <- colnames(meta(corpus))
    varBox1 <- variableListBox(top, vars,
                               title=.gettext("Row variable:"),
                               initialSelection=0)

    varBox2 <- variableListBox(top, vars,
                               title=.gettext("Column variable:"),
                               initialSelection=1)

    radioButtons(name="what",
                 buttons=c("row", "col", "absolute"),
                 labels=c(.gettext("Row %"),
                          .gettext("Column %"),
                          .gettext("Absolute counts")),
                 title=.gettext("Measure:"),
                 right.buttons=FALSE)

    plotFrame <- tkframe(top)

    tclPlotVar <- tclVar(1)
    plotButton <- tkcheckbutton(plotFrame, text=.gettext("Draw plot"), variable=tclPlotVar)

    tclVertVar <- tclVar(0)
    vertButton <- tkcheckbutton(plotFrame, text=.gettext("Vertical bars"), variable=tclVertVar)

    tclStackVar <- tclVar(0)
    stackButton <- tkcheckbutton(plotFrame, text=.gettext("Stacked bars"), variable=tclStackVar)

    tclTitle <- tclVar(.gettext("Distribution of documents by %V1 and %V2"))
    titleEntry <- ttkentry(plotFrame, width=40, textvariable=tclTitle)

    onOK <- function() {
        var1 <- getSelection(varBox1)
        var2 <- getSelection(varBox2)
        what <- tclvalue(whatVariable)
        plot <- tclvalue(tclPlotVar) == 1
        vert <- tclvalue(tclVertVar) == 1
        stack <- tclvalue(tclStackVar) == 1
        title <- tclvalue(tclTitle)

        closeDialog()

        title <- gsub("%V2", tolower(var2), gsub("%V1", tolower(var1), title))

        doItAndPrint(sprintf('absVarFreqs <- table(meta(corpus, c("%s", "%s")))', var1, var2))

        if(what == "row") {
            doItAndPrint("varFreqs <- prop.table(absVarFreqs, 1) * 100")
            ylab <- .gettext("% of documents")
        }
        else if (what == "col") {
            doItAndPrint("varFreqs <- prop.table(absVarFreqs, 2) * 100")
            ylab <- .gettext("% of documents")
        }
        else {
            doItAndPrint("varFreqs <- absVarFreqs")
            ylab <- .gettext("Number of documents")
        }

        # An empty level leads to NAs when computing %
        if(stack && any(is.na(varFreqs))) {
            .Message(.gettext("Cannot plot stacked bars when the table has a null margin."), "error", parent=top)
            stack <- FALSE
        }

        if(plot) {
            if(vert) {
                if(what == "col")
                    doItAndPrint(sprintf('barchart(t(varFreqs), stack=%s, horizontal=FALSE, scales=list(rot=90), ylab="%s", main="%s", auto.key=TRUE)',
                                         stack, ylab, title))
                 else
                    doItAndPrint(sprintf('barchart(varFreqs, stack=%s, horizontal=FALSE, scales=list(rot=90), ylab="%s", main="%s", auto.key=TRUE)',
                                         stack, ylab, title))
            }
            else {
                if(what == "col")
                    doItAndPrint(sprintf('barchart(t(varFreqs), stack=%s, xlab="%s", main="%s", auto.key=TRUE)',
                                         stack, ylab, title))
                 else
                    doItAndPrint(sprintf('barchart(varFreqs, stack=%s, xlab="%s", main="%s", auto.key=TRUE)',
                                         stack, ylab, title))
            }
        }

        if(what == "row")
            doItAndPrint("varFreqs <- addmargins(varFreqs, 2)")
        else if(what == "col")
            doItAndPrint("varFreqs <- addmargins(varFreqs, 1)")
        else
            doItAndPrint("varFreqs <- addmargins(varFreqs)")

        if(what == "row")
            setLastTable("varFreqs", paste(title, .gettext("(row %)")))
        else if(what == "col")
            setLastTable("varFreqs", paste(title, .gettext("(column %)")))
        else
            setLastTable("varFreqs", title)

        doItAndPrint("varFreqs")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="varCrossTableDlg")
    tkgrid(getFrame(varBox1), getFrame(varBox2), sticky="w", pady=6)
    tkgrid(whatFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(.titleLabel(plotFrame, text=.gettext("Plot:")),
           sticky="w", columnspan=2)
    tkgrid(plotButton, sticky="w", columnspan=2)
    tkgrid(vertButton, sticky="w", columnspan=2)
    tkgrid(stackButton, sticky="w", columnspan=2)
    tkgrid(labelRcmdr(plotFrame, text=.gettext("Title:")), titleEntry, sticky="w")
    tkgrid(plotFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=varBox1$listbox)
}

