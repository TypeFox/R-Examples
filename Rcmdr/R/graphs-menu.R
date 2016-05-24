# Graphs menu dialogs

# last modified 2016-02-17 by J. Fox
#  applied patch to improve window behaviour supplied by Milan Bouchet-Valat 2011-09-22

# the following functions improved by Miroslav Ristic 2013-07: barGraph, indexPlot, boxPlot, 
#    DensityPlot, Histogram, pieChart, PlotMeans, QQPlot, scatterPlot, scatterPlotMatrix, Stripchart,
#    Xyplot

indexPlot <- function () {
    defaults <- list(initial.x = NULL, initial.type = "spikes", initial.identify = "auto",
        initial.id.n="2", initial.tab=0, 
        initial.ylab=gettextRcmdr("<auto>"), initial.main=gettextRcmdr("<auto>"))
    dialog.values <- getDialog("indexPlot", defaults)
    initializeDialog(title = gettextRcmdr("Index Plot"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    optionsFrame <- tkframe(optionsTab)    
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Options"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    typeVariable <- tclVar(dialog.values$initial.type)
    styleFrame <- tkframe(optFrame)
    radioButtons(styleFrame, name = "type", buttons = c("spikes", "points"),
        labels = gettextRcmdr(c("Spikes", "Points")), title = gettextRcmdr("Style of plot"),
        initialValue = dialog.values$initial.type)
    identifyPointsFrame <- tkframe(optFrame)
    radioButtons(identifyPointsFrame, name = "identify", buttons = c("auto", "mouse", 
        "not"), labels = gettextRcmdr(c("Automatically", 
            "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"), 
        initialValue = dialog.values$initial.identify)    
    id.n.Var <- tclVar(dialog.values$initial.id.n) 
    npointsSpinner <- tkspinbox(identifyPointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        identify <- tclvalue(identifyVariable)
        id.n <- tclvalue(id.n.Var)
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = indexPlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("indexPlot", list(initial.x = x, initial.type = tclvalue(typeVariable), initial.identify = identify,
            initial.id.n = id.n, initial.tab=tab,
            initial.ylab = tclvalue(ylabVar), initial.main = tclvalue(mainVar)))
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = indexPlot, message = gettextRcmdr("You must select a variable"))
            return()
        }
        type <- if (tclvalue(typeVariable) == "spikes") "h" else "p"
        method <- if (identify == "mouse") "identify" else "y"
        id.n.use <- if (identify == "not") 0 else id.n
        .activeDataSet <- ActiveDataSet()
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
                gettextRcmdr(if (MacOSXP()) 
                    "esc key to exit."
                    else "right button to exit."), sep = ""), icon = "info", 
                type = "ok")
        }
        command <- paste("with(", .activeDataSet, ", indexplot(", x, ", type='", type,
            "', id.method='", method, "', id.n=", id.n.use, ", labels=rownames(", .activeDataSet, ")",
            ylab, main, "))", sep="") # Modification
        if (identify == "mouse") command <- suppressMarkdown(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "indexplot", reset = "indexPlot", apply="indexPlot")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(typeFrame, sticky = "w")
    tkgrid(styleFrame, sticky = "w")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(identifyPointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(identifyPointsFrame, sticky="w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

Histogram <- function () {
    defaults <- list(initial.x = NULL, initial.scale = "frequency",
        initial.bins = gettextRcmdr ("<auto>"), initial.tab=0,
        initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"),
        initial.main=gettextRcmdr("<auto>"), initial.group = NULL)
    dialog.values <- getDialog("Histogram", defaults)
    initializeDialog(title = gettextRcmdr("Histogram"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = Histogram, message = gettextRcmdr("You must select a variable"))
            return()
        }
        bins <- tclvalue(binsVariable)
        opts <- options(warn = -1)
        binstext <- if (bins == gettextRcmdr("<auto>"))
            "\"Sturges\""
        else as.numeric(bins)
        options(opts)
        scale <- tclvalue(scaleVariable)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("Histogram", list (initial.x = x, initial.bins = bins, initial.scale = scale,
            initial.tab=tab, initial.xlab=tclvalue(xlabVar), initial.ylab = tclvalue(ylabVar),
            initial.main = tclvalue(mainVar), initial.group=if (.groups == FALSE) NULL else .groups))
        if (is.null(.groups) || .groups == FALSE) {
            command <- paste("with(", ActiveDataSet(), ", Hist(", x, ', scale="', scale, '", breaks=',
                binstext, ', col="darkgray"', xlab, ylab, main, "))", sep="")
        }
        else{
            command <- paste("with(", ActiveDataSet(), ", Hist(", x, ", groups=", .groups, ', scale="',
                scale, '", breaks=', binstext, ', col="darkgray"', xlab, ylab, main, "))", sep="")
        }
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(Histogram, initialGroup=initial.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups")
        else paste(gettextRcmdr("Plot by:"), initial.group), window=dataTab)
    OKCancelHelp(helpSubject = "Hist", reset = "Histogram", apply="Histogram")
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Options"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    axisFrame <- tkframe(optFrame)
    radioButtons(axisFrame, name = "scale", buttons = c("frequency", "percent",
        "density"), labels = gettextRcmdr(c("Frequency counts",
            "Percentages", "Densities")), title = gettextRcmdr("Axis Scaling"),
        initialValue = dialog.values$initial.scale)
    binsFrame <- tkframe(optFrame)
    binsVariable <- tclVar(dialog.values$initial.bins)
    binsField <- ttkentry(binsFrame, width = "8", textvariable = binsVariable)
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(labelRcmdr(binsFrame, text = gettextRcmdr("Number of bins: ")),
        binsField, sticky = "w")
    tkgrid(binsFrame, sticky = "w")
    tkgrid(scaleFrame, sticky = "w")
    tkgrid(axisFrame, sticky = "w")
    tkgrid.configure(binsField, sticky = "e")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

stemAndLeaf <- function () {
    Library("aplpack")
    defaults <- list(initial.x = NULL, initial.leafs.auto="1", initial.unit = 0,  initial.m = "auto", 
        initial.trim = 1, initial.depths = 1, initial.reverse = 1, initial.style = "Tukey",
        initial.tab=0, initial.group = NULL) 
    dialog.values <- getDialog("stemAndLeaf", defaults)
    initializeDialog(title = gettextRcmdr("Stem and Leaf Display"), 
        #                     preventCrisp = TRUE, use.tabs=TRUE)
        use.tabs=TRUE)
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    displayDigits <- tclVar(formatC(10^dialog.values$initial.unit))
    leafsDigitValue <- tclVar(dialog.values$initial.unit)
    leafsDigitOnce <- TRUE
    onDigits <- function(...) {
        if(leafsDigitOnce) {
            leafsDigitOnce <<- FALSE
            return()
        }
        
        tclvalue(displayDigits) <- formatC(10^as.numeric(tclvalue(leafsDigitValue)), 
            format = "fg", big.mark = ",")
        tclvalue(leafsAutoVariable) <- "0"
    }
    radioButtons(optionsTab, name = "parts", buttons = c("auto", "one", "two", 
        "five"), values = c("auto", "1", "2", "5"), labels = c(gettextRcmdr("Automatic"), 
            "   1", "   2", "   5"), title = gettextRcmdr("Parts Per Stem"), 
        initialValue = dialog.values$initial.m)
    radioButtons(optionsTab, name = "style", buttons = c("Tukey", "bare"), 
        labels = gettextRcmdr(c("Tukey", "Repeated stem digits")), 
        title = gettextRcmdr("Style of Divided Stems"), 
        initialValue = dialog.values$initial.style)
    checkBoxes(optionsTab, frame="otherOptionsFrame", 
        boxes = c("trimOutliers", "showDepths", "reverseNegative"), 
        initialValues = c(dialog.values$initial.trim, dialog.values$initial.depths, dialog.values$initial.reverse),
        labels = gettextRcmdr(c("Trim outliers", "Show depths", "Reverse negative leaves")),
        title=gettextRcmdr("Other Options"))
    leafsFrame <- tkframe(optionsTab)
    leafsDigitValue <- tclVar(dialog.values$initial.unit)
    leafsDigitSlider <- tkscale(leafsFrame, from = -6, to = 6, 
        showvalue = FALSE, variable = leafsDigitValue, resolution = 1, 
        orient = "horizontal", command = onDigits)
    leafsDigitShow <- labelRcmdr(leafsFrame, textvariable = displayDigits, 
        width = 8, justify = "right")
    leafsAutoVariable <- tclVar(dialog.values$initial.leafs.auto)
    leafsDigitCheckBox <- ttkcheckbutton(leafsFrame, text = gettextRcmdr("Automatic"),
        variable = leafsAutoVariable)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        m <- tclvalue(partsVariable)
        style <- tclvalue (styleVariable)
        trim <- tclvalue (trimOutliersVariable)
        depths <- tclvalue (showDepthsVariable)
        reverse <- tclvalue (reverseNegativeVariable)
        unit <- if (tclvalue(leafsAutoVariable) == "1") 
            ""
        else paste(", unit=", 10^as.numeric(tclvalue(leafsDigitValue)), 
            sep = "")
        putDialog ("stemAndLeaf", list(initial.x = x, initial.leafs.auto=tclvalue(leafsAutoVariable),
            initial.unit = as.numeric(tclvalue(leafsDigitValue)),  initial.m = m, 
            initial.trim = trim, initial.depths = depths, initial.reverse = reverse, 
            initial.style = style, initial.tab=tab,
            initial.group=if (.groups == FALSE) NULL else .groups))
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = stemAndLeaf, message = gettextRcmdr("You must select a variable"))
            return()
        }
        trim <- if (tclvalue(trimOutliersVariable) == "1") 
            ""
        else ", trim.outliers=FALSE"
        depths <- if (tclvalue(showDepthsVariable) == "1") 
            ""
        else ", depths=FALSE"
        reverse <- if (tclvalue(reverseNegativeVariable) == "1") 
            ""
        else ", reverse.negative.leaves=FALSE"
        m <- if (tclvalue(partsVariable) == "auto") 
            ""
        else paste(", m=", tclvalue(partsVariable), sep = "")
        style <- if (tclvalue(styleVariable) == "Tukey") 
            ""
        else ", style=\"bare\""
        command <- if (is.null(.groups) || .groups == FALSE){
            paste("with(", ActiveDataSet(), ", stem.leaf(",  
                x, style, unit, m, trim, depths, reverse, ", na.rm=TRUE))", 
                sep = "")
        }
        else {
            levels <- levels(eval(parse(text=paste(ActiveDataSet(), "$", .groups, sep="")), 
                envir=.GlobalEnv))
            paste("with(", ActiveDataSet(),
                ", stem.leaf.backback(", x,"[", .groups, ' == "', levels[1], '"], ',
                x,"[", .groups, ' == "', levels[2], '"]',
                style, unit, m, trim, depths, reverse, ", na.rm=TRUE))", 
                sep = "")
        }
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    groupsBox(stemAndLeaf, variables=TwoLevelFactors(), initialGroup=initial.group,
        label=gettextRcmdr("Plot back-to-back by:"),
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot back-to-back by")
        else paste(gettextRcmdr("Plot back-to-back by:"), initial.group), 
        errorText=gettextRcmdr("There are no two-level factors in the active data set."), 
        window=dataTab)
    OKCancelHelp(helpSubject = "stem.leaf", reset = "stemAndLeaf", apply = "stemAndLeaf")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(labelRcmdr(leafsFrame, text = gettextRcmdr("Leafs Digit:  "), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), leafsDigitCheckBox, 
        labelRcmdr(leafsFrame, text = gettextRcmdr("  or set:"), 
            fg = "red"), leafsDigitShow, leafsDigitSlider, sticky = "w")
    tkgrid(partsFrame, sticky = "w")
    tkgrid(styleFrame, sticky = "w")
    tkgrid(otherOptionsFrame, sticky="w")
    tkgrid(leafsFrame, sticky = "w")
    tkgrid(groupsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

boxPlot <- function () {
    defaults <- list(initial.x = NULL, initial.identify = "y", initial.group=NULL,
        initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"), 
        initial.main=gettextRcmdr("<auto>"), initial.tab=0)
    dialog.values <- getDialog("boxPlot", defaults)
    initializeDialog(title = gettextRcmdr("Boxplot"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Identify Outliers"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    styleFrame <- tkframe(optFrame)
    radioButtons(optFrame, name = "identify", buttons = c("y", "identify", "none"),
        labels = gettextRcmdr(c("Automatically", "With mouse", "No")),
        initialValue = dialog.values$initial.identify, )
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        identifyPoints <- tclvalue(identifyVariable)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("boxPlot", list(initial.x = x, initial.identify = identifyPoints,
            initial.group=if (.groups == FALSE) NULL else .groups,
            initial.xlab=tclvalue(xlabVar), initial.ylab=tclvalue(ylabVar),
            initial.main=tclvalue(mainVar), initial.tab=tab))
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = boxPlot, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        var <- paste(.activeDataSet, "$", x, sep = "")
        if (identifyPoints == "identify")
            RcmdrTkmessageBox(title = "Identify Points",
                message = paste(gettextRcmdr("Use left mouse button to identify points,\n"),
                    gettextRcmdr(if (MacOSXP()) "esc key to exit."
                        else "right button to exit."), sep = ""),
                icon = "info", type = "ok")
        if (is.null(.groups) || .groups == FALSE) {
            command <- paste("Boxplot( ~ ", x, ", data=", .activeDataSet, ', id.method="',
                identifyPoints, '"', ylab, main, ')', sep="")
            if (identifyPoints == "identify") command <- suppressMarkdown(command)
            
            doItAndPrint(command)
        }
        else {
            command <- paste("Boxplot(", x, "~", .groups, ", data=", .activeDataSet,
                ', id.method="', identifyPoints, '"', xlab, ylab, main, ')', sep = "")
            if (identifyPoints == "identify") command <- suppressMarkdown(command)
            doItAndPrint(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(boxPlot, initialGroup=initial.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups")
        else paste(gettextRcmdr("Plot by:"), initial.group), window=dataTab)
    OKCancelHelp(helpSubject = "Boxplot", reset = "boxPlot", apply="boxPlot")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(identifyFrame, stick = "w")
    tkgrid(styleFrame, stick="w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

DotPlot <- function () {
    defaults <- list(initial.x = NULL, initial.use.bins = 0,
                     initial.bins = gettextRcmdr ("<auto>"), initial.tab=0,
                     initial.xlab=gettextRcmdr("<auto>"),
                     initial.group = NULL)
    dialog.values <- getDialog("DotPlot", defaults)
    initializeDialog(title = gettextRcmdr("Dot Plot"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"),
                            initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = DotPlot, message = gettextRcmdr("You must select a variable"))
            return()
        }
        use.bins <- tclvalue(useBinsVariable) == "1"
        bins <- tclvalue(binsVariable)
        opts <- options(warn = -1)
        binstext <- if (bins == gettextRcmdr("<auto>"))
            "\"Sturges\""
        else as.numeric(bins)
        binstext <- if (use.bins) paste(', breaks=', binstext, sep="") else ""
        options(opts)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        putDialog ("DotPlot", list (initial.x = x, initial.bins = bins, initial.use.bins = if (use.bins) "1" else "0",
                                    initial.tab=tab, initial.xlab=tclvalue(xlabVar), 
                                    initial.group=if (.groups == FALSE) NULL else .groups))
        if (is.null(.groups) || .groups == FALSE) {
            command <- paste("with(", ActiveDataSet(), ", Dotplot(", x, ', bin=', use.bins,
                             binstext, xlab, "))", sep="")
        }
        else{
            command <- paste("with(", ActiveDataSet(), ", Dotplot(", x, ", by=", .groups, ', bin=',
                             use.bins, binstext, xlab, "))", sep="")
        }
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(DotPlot, initialGroup=initial.group,
              initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups")
              else paste(gettextRcmdr("Plot by:"), initial.group), window=dataTab)
    OKCancelHelp(helpSubject = "Dotplot", reset = "DotPlot", apply="DotPlot", helpPackage="RcmdrMisc")
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Options"),
                                                                font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame)
    xlabVar <- tclVar(dialog.values$initial.xlab)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    binsFrame <- tkframe(optFrame)
    useBinsVariable <- tclVar(dialog.values$initial.use.bins)
    useBinsCheckBox <- ttkcheckbutton(binsFrame, text = gettextRcmdr("Bin variable"),
                                      variable = useBinsVariable)
    binsVariable <- tclVar(dialog.values$initial.bins)
    binsField <- ttkentry(binsFrame, width = "8", textvariable = binsVariable)
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(useBinsCheckBox, sticky="w")
    tkgrid(labelRcmdr(binsFrame, text = gettextRcmdr("Number of bins: ")),
           binsField, sticky = "w")
    tkgrid(binsFrame, sticky = "w")
    tkgrid.configure(binsField, sticky = "e")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

scatterPlot <- function () {
    defaults <- list(initial.x = NULL, initial.y = NULL, initial.jitterx = 0, initial.jittery = 0,
        initial.logstringx = 0, initial.logstringy = 0, initial.box = 0,
        initial.line = 0, initial.smooth = 0, initial.spread = 0, initial.span = 50,
        initial.ellipse=0, initial.levels=".5, .9",
        initial.subset = gettextRcmdr ("<all valid cases>"), initial.ylab = gettextRcmdr ("<auto>"),
        initial.xlab = gettextRcmdr("<auto>"), initial.pch = gettextRcmdr("<auto>"),
        initial.cexValue = 1, initial.cex.axisValue = 1, initial.cex.labValue = 1, initialGroup=NULL, initial.lines.by.group=1,
        initial.identify=gettextRcmdr("not"), initial.identify.points="2", initial.tab=0,
        initial.main=gettextRcmdr("<auto>"), initial.legend.pos="above")
    dialog.values <- getDialog("scatterPlot", defaults)
    initial.group <- dialog.values$initial.group
    .linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    initializeDialog(title = gettextRcmdr("Scatterplot"), use.tabs=TRUE)
    .numeric <- Numeric()
    xBox <- variableListBox(dataTab, .numeric, title = gettextRcmdr("x-variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    yBox <- variableListBox(dataTab, .numeric, title = gettextRcmdr("y-variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.y, "numeric"))
    optionsParFrame <- tkframe(optionsTab)
    parFrame <- ttklabelframe(optionsParFrame, labelwidget=tklabel(optionsParFrame, text=gettextRcmdr("Plot Labels and Points"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    checkBoxes(window = optionsParFrame, frame = "optionsFrame",
        boxes = c("jitterX", "jitterY", "logX", "logY",
            "boxplots", "lsLine", "smoothLine", "spread"), initialValues = c(
                dialog.values$initial.jitterx, dialog.values$initial.jittery,
                dialog.values$initial.logstringx, dialog.values$initial.logstringy,
                dialog.values$initial.box, dialog.values$initial.line, dialog.values$initial.smooth,
                dialog.values$initial.spread),labels = gettextRcmdr(c(
                    "Jitter x-variable", "Jitter y-variable", "Log x-axis",
                    "Log y-axis", "Marginal boxplots", "Least-squares line",
                    "Smooth line", "Show spread")), title = gettextRcmdr("Plot Options"), ttk=TRUE)
    sliderValue <- tclVar(dialog.values$initial.span)
    sliderFrame <- tkframe(optionsFrame)
    slider <- tkscale(sliderFrame, from = 5, to = 100, showvalue = TRUE,
        variable = sliderValue, resolution = 5, orient = "horizontal")
    ellipseVariable <- tclVar(dialog.values$initial.ellipse)
    ellipseFrame <- tkframe(optionsFrame)
    ellipseCheckBox <- ttkcheckbutton(ellipseFrame, variable=ellipseVariable)
    levelsVar <- tclVar(dialog.values$initial.levels)
    levelsFrame <- tkframe(optionsFrame)
    levelsEntry <- ttkentry(levelsFrame, width = "10", textvariable = levelsVar)
    radioButtons(window=optionsFrame, name = "identify", buttons = c("auto", "mouse",
        "not"), labels = gettextRcmdr(c("Automatically",
            "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"),
        initialValue = dialog.values$initial.identify)
    id.n.Var <- tclVar(dialog.values$initial.identify.points)
    npointsFrame <- tkframe(optionsFrame)
    npointsSpinner <- tkspinbox(npointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)
    subsetBox(dataTab, subset.expression = dialog.values$initial.subset)
    tkbind(subsetEntry, "<FocusIn>", function() tkselection.clear(subsetEntry))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkbind(xlabEntry, "<FocusIn>", function() tkselection.clear(xlabEntry))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    pchVar <- tclVar(dialog.values$initial.pch)
    pchEntry <- ttkentry(parFrame, width = 25, textvariable = pchVar)
    cexValue <- tclVar(dialog.values$initial.cexValue)
    cex.axisValue <- tclVar(dialog.values$initial.cex.axisValue)
    cex.labValue <- tclVar(dialog.values$initial.cex.labValue)
    cexSlider <- tkscale(parFrame, from = 0.5, to = 2.5, showvalue = TRUE,
        variable = cexValue, resolution = 0.1, orient = "horizontal")
    cex.axisSlider <- tkscale(parFrame, from = 0.5, to = 2.5,
        showvalue = TRUE, variable = cex.axisValue, resolution = 0.1,
        orient = "horizontal")
    cex.labSlider <- tkscale(parFrame, from = 0.5, to = 2.5,
        showvalue = TRUE, variable = cex.labValue, resolution = 0.1,
        orient = "horizontal")
    radioButtons(window=parFrame, name="legendPosition", buttons=c("above", "topleft", "topright",
        "bottomleft", "bottomright"), labels=gettextRcmdr(c("Above plot", "Top left", "Top right", "Bottom left",
            "Bottom right")), title=gettextRcmdr("Legend Position"),
        initialValue=dialog.values$initial.legend.pos)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        jitter <- if ("1" == tclvalue(jitterXVariable) && "1" ==
                tclvalue(jitterYVariable))
            ", jitter=list(x=1, y=1)"
        else if ("1" == tclvalue(jitterXVariable))
            ", jitter=list(x=1)"
        else if ("1" == tclvalue(jitterYVariable))
            ", jitter=list(y=1)"
        else ""
        logstring <- ""
        if ("1" == tclvalue(logXVariable))
            logstring <- paste(logstring, "x", sep = "")
        if ("1" == tclvalue(logYVariable))
            logstring <- paste(logstring, "y", sep = "")
        identify <- tclvalue(identifyVariable)
        id.n <- tclvalue(id.n.Var)
        identify.text <- switch(identify,
            auto = paste(", id.method='mahal', id.n =", id.n),
            mouse = ", id.method='identify'",
            not = "")
        box <- tclvalue(boxplotsVariable)
        line <- tclvalue(lsLineVariable)
        smooth <-  tclvalue(smoothLineVariable)
        spread <- tclvalue(spreadVariable)
        span <- as.numeric(tclvalue(sliderValue))
        initial.ellipse <- tclvalue(ellipseVariable)
        ellipse <- as.character(initial.ellipse == "1")
        save.levels <- levels <- tclvalue(levelsVar)
        levels <- gsub(",", " ", levels)
        levels <- paste("c(", gsub('[ ]+', ", ", levels), ")", sep="")
        res <- try(is.numeric(eval(parse(text=levels))), silent=TRUE)
        if (class(res) == "try-error" || !res){
            errorCondition(recall = scatterPlot, message = 
                    gettextRcmdr("Levels for ellipses must be numeric values\n  separated by spaces or commas."))
            return()
        }
        initial.subset <- subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>"))
            ""
        else paste(", subset=", subset, sep = "")
        cex.axis <- as.numeric(tclvalue(cex.axisValue))
        cex <- as.numeric(tclvalue(cexValue))
        cex.lab <- as.numeric(tclvalue(cex.labValue))
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        pchVal <- gsub(" ", ",", tclvalue(pchVar))
        legend.pos <- tclvalue(legendPositionVariable)
        closeDialog()
        if ("" == pchVal) {
            errorCondition(recall = scatterPlot, message = gettextRcmdr("No plotting characters."))
            return()
        }
        pch <- if (trim.blanks(pchVal) == gettextRcmdr("<auto>"))
            ""
        else paste(", pch=c(", pchVal, ")", sep = "")
        if (length(x) == 0 || length(y) == 0) {
            errorCondition(recall = scatterPlot, message = gettextRcmdr("You must select two variables"))
            return()
        }
        if (x == y) {
            errorCondition(recall = scatterPlot, message = gettextRcmdr("x and y variables must be different"))
            return()
        }
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = scatterPlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        putDialog ("scatterPlot", list (initial.x = x, initial.y = y, initial.jitterx = tclvalue(jitterXVariable),
            initial.jittery = tclvalue(jitterYVariable), initial.logstringx = tclvalue(logXVariable),
            initial.logstringy = tclvalue(logYVariable),  initial.box = box,
            initial.line = line, initial.smooth = smooth, initial.spread = spread,
            initial.span = span, initial.ellipse = initial.ellipse, initial.levels = save.levels,
            initial.subset = initial.subset, initial.xlab = tclvalue(xlabVar),
            initial.ylab = tclvalue(ylabVar), initial.cexValue = tclvalue(cexValue),
            initial.cex.axisValue = tclvalue(cex.axisValue), initial.cex.labValue = tclvalue(cex.labValue),
            initial.pch = pchVal, initial.group=if (.groups == FALSE) NULL else .groups,
            initial.lines.by.group=if (.linesByGroup) 1 else 0, initial.identify=identify, initial.identify.points=id.n,
            initial.tab=tab, initial.main=tclvalue(mainVar), initial.legend.pos=legend.pos)
        )
        .activeDataSet <- ActiveDataSet()
        log <- if (logstring != "")
            paste(", log=\"", logstring, "\"", sep = "")
        else ""
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = gettextRcmdr("Identify Points"), message = paste(gettextRcmdr("Use left mouse button to identify points,\n"),
                gettextRcmdr(if (MacOSXP())
                    "esc key to exit."
                    else "right button to exit."), sep = ""), icon = "info",
                type = "ok")
        }
        box <- if ("1" == tclvalue(boxplotsVariable))
            "'xy'"
        else "FALSE"
        line <- if ("1" == tclvalue(lsLineVariable))
            "lm"
        else "FALSE"
        smooth <- as.character("1" == tclvalue(smoothLineVariable))
        spread <- as.character("1" == tclvalue(spreadVariable))
        cex <- if (cex == 1)
            ""
        else paste(", cex=", cex, sep = "")
        cex.axis <- if (cex.axis == 1)
            ""
        else paste(", cex.axis=", cex.axis, sep = "")
        cex.lab <- if (cex.lab == 1)
            ""
        else paste(", cex.lab=", cex.lab, sep = "")
        if (.groups == FALSE) {
            command <- paste("scatterplot(", y, "~", x, log,
                ", reg.line=", line, ", smooth=", smooth, ", spread=",
                spread, identify.text, ", boxplots=", box, ", span=", span/100, 
                ", ellipse=", ellipse, ", levels=", levels,
                jitter, xlab, ylab, main, cex, cex.axis,
                cex.lab, pch, ", data=", .activeDataSet, subset,
                ")", sep = "")
            if (identify == "mouse") command <- suppressMarkdown(command)
            doItAndPrint(command)
        }
        else {
            command <- paste("scatterplot(", y, "~", x, " | ",
                .groups, log, ", reg.line=", line, ", smooth=", smooth,
                ", spread=", spread, identify.text, ", boxplots=", box,
                ", span=", span/100, 
                ", ellipse=", ellipse, ", levels=", levels,
                jitter, xlab, ylab, main, cex,
                cex.axis, cex.lab, pch, ", by.groups=", .linesByGroup,
                if (legend.pos != "above") paste0(', legend.coords="', legend.pos, '"'),
                ", data=", .activeDataSet, subset, ")", sep = "")
            if (identify == "mouse") command <- suppressMarkdown(command)
            doItAndPrint(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group), window=dataTab)
    OKCancelHelp(helpSubject = "scatterplot", reset = "scatterPlot", apply="scatterPlot")
    tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw", padx=6, pady=c(6, 0))
    tkgrid(labelRcmdr(sliderFrame, text = gettextRcmdr("Span for smooth")),
        slider, sticky = "swe", padx=6, pady=6)
    tkgrid(sliderFrame, sticky="w")
    tkgrid(ellipseCheckBox, labelRcmdr(ellipseFrame, text=gettextRcmdr("Plot concentration ellipse(s)")), sticky="w")
    tkgrid(ellipseFrame, sticky="w", pady=6)
    tkgrid(labelRcmdr(levelsFrame, text=gettextRcmdr("Concentration levels: ")), levelsEntry, sticky="w")
    tkgrid(levelsFrame, sticky="w")
    tkgrid(identifyFrame, sticky="w", pady=6)
    tkgrid(labelRcmdr(npointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(npointsFrame, sticky="w")
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Plotting characters")),
        pchEntry, stick = "we", padx=6, pady=6)
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Point size")),
        cexSlider, sticky = "we", padx=6, pady=6)
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis text size")),
        cex.axisSlider, sticky = "we", padx=6, pady=6)
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis-labels text size")),
        cex.labSlider, sticky = "we", padx=6, pady=6)
    tkgrid(legendPositionFrame, stick="w", pady=6)
    tkgrid(optionsFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsParFrame, sticky = "we")
    tkgrid(ttklabel(dataTab, text=""))
    tkgrid(groupsFrame, sticky = "we", padx=6)
    tkgrid(ttklabel(dataTab, text=""))
    tkgrid(subsetFrame, sticky = "we", padx=6, pady=c(0, 6))
    tkgrid(labelRcmdr(top, text = " "), padx=6)
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

scatterPlotMatrix <- function () {
    defaults <- list(initial.variables = NULL, initial.line = 0, initial.smooth = 0, initial.spread = 0,
        initial.span = 50, initial.ellipse=0, initial.levels=".5, .9",
        initial.diag = "density", initial.subset = gettextRcmdr ("<all valid cases>"),
        initialGroup=NULL, initial.lines.by.group=1, initial.id.n="0", initial.tab=0,
        initial.main=gettextRcmdr("<auto>"))
    dialog.values <- getDialog("scatterPlotMatrix", defaults)
    initial.group <- dialog.values$initial.group
    .linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    initializeDialog(title = gettextRcmdr("Scatterplot Matrix"), use.tabs=TRUE)
    variablesBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Select variables (three or more)"),
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- tkframe(optionsFrame) #ttklabelframe(optionsFrame, text = gettextRcmdr("Options"))
    parFrame <- tkframe(optionsFrame) #ttklabelframe(optionsFrame, text = gettextRcmdr("Plotting Parameters"))
    checkBoxes(window = optFrame, frame = "otherFrame", boxes = c("lsLine", "smoothLine", "spread"),
        initialValues = c(dialog.values$initial.line, dialog.values$initial.smooth,
            dialog.values$initial.spread), labels = gettextRcmdr(c("Least-squares lines",
                "Smooth lines", "Show spread")),
        title=gettextRcmdr("Other Options"))
    identifyFrame <- tkframe(optFrame)
    id.n.Var <- tclVar(dialog.values$initial.id.n)
    npointsSpinner <- tkspinbox(identifyFrame, from=0, to=10, width=2, textvariable=id.n.Var)
    sliderValue <- tclVar(dialog.values$initial.span)
    sliderFrame <- tkframe(optFrame)
    slider <- tkscale(sliderFrame, from = 5, to = 100, showvalue = TRUE,
        variable = sliderValue, resolution = 5, orient = "horizontal")
    ellipseVariable <- tclVar(dialog.values$initial.ellipse)
    ellipseFrame <- tkframe(optionsFrame)
    ellipseCheckBox <- ttkcheckbutton(ellipseFrame, variable=ellipseVariable)
    levelsVar <- tclVar(dialog.values$initial.levels)
    levelsFrame <- tkframe(optionsFrame)
    levelsEntry <- ttkentry(levelsFrame, width = "10", textvariable = levelsVar)
    radioButtons(window = optFrame, name = "diagonal", buttons = c("density", "histogram",
        "boxplot", "oned", "qqplot", "none"), labels = gettextRcmdr(c("Density plots",
            "Histograms", "Boxplots", "One-dimensional scatterplots",
            "Normal QQ plots", "Nothing (empty)")), title = gettextRcmdr("On Diagonal"),
        initialValue = dialog.values$initial.diag)
    subsetBox(dataTab, subset.expression = dialog.values$initial.subset)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 3) {
            errorCondition(recall = scatterPlotMatrix, message = gettextRcmdr("Fewer than 3 variable selected."))
            return()
        }
        line <- if ("1" == tclvalue(lsLineVariable))
            "lm"
        else "FALSE"
        smooth <- as.character("1" == tclvalue(smoothLineVariable))
        spread <- as.character("1" == tclvalue(spreadVariable))
        id.n <- tclvalue(id.n.Var)
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = scatterPlotMatrix,
                message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        span <- as.numeric(tclvalue(sliderValue))
        initial.ellipse <- tclvalue(ellipseVariable)
        ellipse <- as.character(initial.ellipse == "1")
        save.levels <- levels <- tclvalue(levelsVar)
        levels <- gsub(",", " ", levels)
        levels <- paste("c(", gsub('[ ]+', ", ", levels), ")", sep="")
        res <- try(is.numeric(eval(parse(text=levels))), silent=TRUE)
        if (class(res) == "try-error" || !res){
            errorCondition(recall = scatterPlotMatrix, message = 
                    gettextRcmdr("Levels for ellipses must be numeric values\n  separated by spaces or commas."))
            return()
        }
        diag <- as.character(tclvalue(diagonalVariable))
        initial.subset <- subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
        else paste(", subset=", subset, sep="")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog("scatterPlotMatrix", list(initial.variables = variables, initial.line = tclvalue (lsLineVariable),
            initial.smooth = tclvalue(smoothLineVariable),initial.spread = tclvalue (spreadVariable),
            initial.span = span, initial.ellipse = initial.ellipse, initial.levels = save.levels,
            initial.diag = diag, initial.subset = initial.subset,
            initial.group=if (.groups == FALSE) NULL else .groups,
            initial.lines.by.group=if (.linesByGroup) 1 else 0,
            initial.id.n=id.n,
            initial.tab=tab, initial.main=tclvalue(mainVar)))
        .activeDataSet <- ActiveDataSet()
        if (.groups == FALSE) {
            command <- paste("scatterplotMatrix(~", paste(variables,
                collapse = "+"), ", reg.line=", line, ", smooth=",
                smooth, ", spread=", spread, ", span=", span/100, 
                ", ellipse=", ellipse, ", levels=", levels,
                ", id.n=", id.n,
                ", diagonal = '", diag, "', data=", .activeDataSet,
                subset, main, ")", sep = "")
            logger(command)
            justDoIt(command)
        }
        else {
            command <- paste("scatterplotMatrix(~", paste(variables,
                collapse = "+"), " | ", .groups, ", reg.line=",
                line, ", smooth=", smooth, ", spread=", spread,
                ", span=", span/100,  
                ", ellipse=", ellipse, ", levels=", levels,
                ", id.n=", id.n,
                ", diagonal= '", diag, "', by.groups=",
                .linesByGroup, ", data=", .activeDataSet, subset, main,
                ")", sep = "")
            logger(command)
            justDoIt(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlotMatrix, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group),
        window=dataTab)
    OKCancelHelp(helpSubject = "scatterplotMatrix", reset = "scatterPlotMatrix", apply = "scatterPlotMatrix")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(diagonalFrame, sticky = "w")
    tkgrid(otherFrame, sticky = "w")
    tkgrid(labelRcmdr(sliderFrame, text = gettextRcmdr("Span for smooth")),
        slider, sticky = "w")
    tkgrid(sliderFrame, sticky="w")
    mainVar <- tclVar(dialog.values$initial.main)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(identifyFrame,
        text = gettextRcmdr("Number of points to identify  \nin each panel and group ")),
        npointsSpinner, sticky="nw")
    tkgrid(identifyFrame, sticky="w", columnspan=2, pady=6)
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(ellipseCheckBox, labelRcmdr(ellipseFrame, text=gettextRcmdr("Plot concentration ellipse(s)")), sticky="w")
    tkgrid(ellipseFrame, sticky="w", padx=6)
    tkgrid(labelRcmdr(levelsFrame, text=gettextRcmdr("Concentration levels: ")), levelsEntry, sticky="w")
    tkgrid(levelsFrame, sticky="w", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(subsetFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

barGraph <- function () {
    defaults <- list (initial.variable = NULL, initial.xlab=gettextRcmdr("<auto>"),
                      initial.ylab=gettextRcmdr("<auto>"), initial.main=gettextRcmdr("<auto>"),
                      initial.group=NULL, initial.style="divided", initial.legend="topright",
                      initial.tab=0)
    dialog.values <- getDialog ("barGraph", defaults)
    initializeDialog(title = gettextRcmdr("Bar Graph"), use.tabs=TRUE)
    optionsFrame <- tkframe(optionsTab)
    optionsFrame2 <- tkframe(optionsTab)
    variablesFrame <- tkframe(dataTab)
    variableBox <- variableListBox(variablesFrame, Factors(), title = gettextRcmdr("Variable (pick one)"),
                                   initialSelection = varPosn (dialog.values$initial.variable, "factor"))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
                                                                font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        variable <- getSelection(variableBox)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            paste(", xlab=\"", variable, "\"", sep = "")
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            paste(", ylab=\"Frequency\"", sep = "")
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        style <- tclvalue(styleVariable)
        legend <- tclvalue(legendVariable)
        putDialog ("barGraph", list(initial.variable = variable, initial.xlab=tclvalue(xlabVar),
                                    initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar), 
                                    initial.group=if (.groups == FALSE) NULL else .groups,
                                    initial.style=style, initial.legend=legend, initial.tab=tab))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = barGraph, message = gettextRcmdr("You must select a variable"))
            return()
        }
        command <- if (.groups == FALSE){
            paste("with(", ActiveDataSet(), ", Barplot(", 
                  variable, xlab, ylab, main, "))",
                  sep = "")
        }
        else {
            if (.groups == variable) {
                errorCondition(recall=barGraph, message=gettextRcmdr("plotting and grouping variables must be different"))
                return()
            }
            paste("with(", ActiveDataSet(), ", Barplot(", 
                   variable, ", by=", .groups, ', style="', style, 
                   '", legend.pos="', legend, '"',
                   xlab, ylab, main, "))",
                   sep = "")
        }
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    initial.group <- dialog.values$initial.group
    groupsBox(barGraph, initialGroup=initial.group, 
              initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") 
              else paste(gettextRcmdr("Plot by:"), initial.group), window=variablesFrame)
    radioButtons(optionsFrame2, name = "style", buttons = c("divided", "parallel"), 
                 labels = gettextRcmdr(c("Divided (stacked)", "Side-by-side (parallel)")), 
                 title = gettextRcmdr("Style of Group Bars"),
                 initialValue = dialog.values$initial.style)
    radioButtons(optionsFrame2, name = "legend", 
                 buttons = c("topright", "top", "topleft"), 
                 labels = gettextRcmdr(c("Right", "Center", "Left")), 
                 title = gettextRcmdr("Position of Legend"),
                 initialValue = dialog.values$initial.legend)
    OKCancelHelp(helpSubject = "Barplot", reset = "barGraph", apply = "barGraph")
    tkgrid(getFrame(variableBox), sticky="w")
    tkgrid(tklabel(variablesFrame, text=""))
    tkgrid(groupsFrame, sticky="w")
    tkgrid(styleFrame, sticky="w")
    tkgrid(legendFrame, sticky="w")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(parFrame, sticky = "nw")
    tkgrid(optionsFrame2, optionsFrame, sticky = "nw")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

pieChart <- function () {
    Library("colorspace")
    defaults <- list (initial.variable = NULL, initial.xlab=gettextRcmdr("<auto>"),
        initial.ylab=gettextRcmdr("<auto>"), initial.main=gettextRcmdr("<auto>"))
    dialog.values <- getDialog ("pieChart", defaults)
    initializeDialog(title = gettextRcmdr("Pie Chart"))
    optionsFrame <- tkframe(top)
    variableBox <- variableListBox(optionsFrame, Factors(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.variable, "factor"))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        variable <- getSelection(variableBox)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            paste(", xlab=\"\"", sep = "")
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            paste(", ylab=\"\"", sep = "")
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            paste(", main=\"", variable, "\"", sep = "")
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("pieChart", list (initial.variable = variable, initial.xlab=tclvalue(xlabVar),
            initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar)))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = pieChart, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        command <- (paste("with(", .activeDataSet, ", pie(table(", 
            variable, "), labels=levels(",
            variable, ")", xlab, ylab, main, ", col=rainbow_hcl(length(levels(",
            variable, ")))))", sep = ""))
        logger(command)
        justDoIt(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "pie", reset = "pieChart", apply = "pieChart")
    tkgrid(getFrame(variableBox), parFrame, sticky = "nw")
    tkgrid(parFrame, sticky = "nw")
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(grid.buttons=TRUE)
}

linePlot <- function () {
    defaults <- list(initial.x = NULL, initial.y = NULL) 
    dialog.values <- getDialog("linePlot", defaults)
    initializeDialog(title = gettextRcmdr("Line Plot"))
    variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("x variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    yBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("y variables (pick one or more)"), 
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.y, "numeric"))
    onOK <- function() {
        y <- getSelection(yBox)
        x <- getSelection(xBox)
        closeDialog()
        if (0 == length(x)) {
            errorCondition(recall = linePlot, message = gettextRcmdr("No x variable selected."))
            return()
        }
        if (0 == length(y)) {
            errorCondition(recall = linePlot, message = gettextRcmdr("No y variables selected."))
            return()
        }
        if (is.element(x, y)) {
            errorCondition(recall = linePlot, message = gettextRcmdr("x and y variables must be different."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        .x <- na.omit(eval(parse(text = paste(.activeDataSet, 
            "$", x, sep = "")), envir = .GlobalEnv))
        if (!identical(order(.x), seq(along.with = .x))) {
            response <- tclvalue(RcmdrTkmessageBox(message = gettextRcmdr("x-values are not in order.\nContinue?"), 
                icon = "warning", type = "okcancel", default = "cancel"))
            if (response == "cancel") {
                onCancel()
                return()
            }
        }
        putDialog ("linePlot", list(initial.x = x, initial.y = y))
        command <- paste("with(", .activeDataSet, ", lineplot(", x, ", ", paste(y, collapse=", "), "))", sep="")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "lineplot", reset = "linePlot", apply = "linePlot")
    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text = "    "), 
        getFrame(yBox), sticky = "nw")
    tkgrid(variablesFrame, sticky = "nw")
    tkgrid(buttonsFrame, stick = "w")
    dialogSuffix()
}

QQPlot <- function () {
    # this function modified by Martin Maechler
    defaults <- list(initial.x = NULL, initial.dist = "norm", initial.df = "",
        initial.chisqdf = "", initial.fdf1 = "", initial.fdf2 = "", initial.othername = "",
        initial.otherparam = "", initial.identify = "auto", initial.id.n="2",
        initial.tab=0, initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"),
        initial.main=gettextRcmdr("<auto>"))
    dialog.values <- getDialog("QQPlot", defaults)
    initializeDialog(title = gettextRcmdr("Quantile-Comparison (QQ) Plot"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Options"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    identifyPointsFrame <- tkframe(optFrame)
    radioButtons(identifyPointsFrame, name = "identify", buttons = c("auto", "mouse",
        "not"), labels = gettextRcmdr(c("Automatically",
            "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"),
        initialValue = dialog.values$initial.identify)
    id.n.Var <- tclVar(dialog.values$initial.id.n)
    npointsSpinner <- tkspinbox(identifyPointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        initial.dist <-dist <- tclvalue(distVariable)
        tdf <- tclvalue(tDfVariable)
        chisqdf <- tclvalue(chisqDfVariable)
        fdf1 <- tclvalue(FDf1Variable)
        fdf2 <- tclvalue(FDf2Variable)
        othername <- tclvalue(otherNameVariable)
        otherparam <- tclvalue(otherParamsVariable)
        id.n <- tclvalue(id.n.Var)
        identify <- tclvalue(identifyVariable)
        method <- if (identify == "mouse") "identify" else "y"
        id.n.use <- if (identify == "not") 0 else id.n
        closeDialog()
        if (0 == length(x)) {
            errorCondition(recall = QQPlot, message = gettextRcmdr("You must select a variable."))
            return()
        }
        save <- options(warn = -1)
        on.exit(save)
        retryMe <- function(msg) {
            Message(message = msg, type = "error")
            QQPlot()
        }
        switch(dist, norm = {
            args <- "dist=\"norm\""
        }, t = {
            df <- tclvalue(tDfVariable)
            df.num <- as.numeric(df)
            if (is.na(df.num) || df.num < 1) {
                retryMe(gettextRcmdr("df for t must be a positive number."))
                return()
            }
            args <- paste("dist=\"t\", df=", df, sep = "")
        }, chisq = {
            df <- tclvalue(chisqDfVariable)
            df.num <- as.numeric(df)
            if (is.na(df.num) || df.num < 1) {
                retryMe(gettextRcmdr("df for chi-square must be a positive number."))
                return()
            }
            args <- paste("dist=\"chisq\", df=", df, sep = "")
        }, f = {
            df1 <- tclvalue(FDf1Variable)
            df2 <- tclvalue(FDf2Variable)
            df.num1 <- as.numeric(df1)
            df.num2 <- as.numeric(df2)
            if (is.na(df.num1) || df.num1 < 1 || is.na(df.num2) ||
                    df.num2 < 1) {
                retryMe(gettextRcmdr("numerator and denominator \ndf for F must be positive numbers."))
                return()
            }
            args <- paste("dist=\"f\", df1=", df1, ", df2=",
                df2, sep = "")
        }, {
            dist <- tclvalue(otherNameVariable)
            params <- tclvalue(otherParamsVariable)
            args <- paste("dist=\"", dist, "\", ", params, sep = "")
        })
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = QQPlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("QQPlot", list (initial.x = x, initial.dist = initial.dist,
            initial.identify = identify, initial.df = tdf, initial.chisqdf = chisqdf,
            initial.fdf1 = fdf1, initial.fdf2 = fdf2, initial.othername = othername,
            initial.otherparam = otherparam, initial.identify = identify, initial.id.n=id.n,
            initial.tab=tab, initial.xlab=tclvalue(xlabVar),
            initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar)))
        .activeDataSet <- ActiveDataSet()
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"),
                gettextRcmdr(if (MacOSXP())
                    "esc key to exit."
                    else "right button to exit."), sep = ""), icon = "info",
                type = "ok")
        }
        command <- paste("with(", .activeDataSet, ", qqPlot", "(", 
            x, ", ", args, ', id.method="', method, '", id.n=', id.n.use, ", labels=rownames(", .activeDataSet, ")", xlab,
            ylab, main, "))", sep = "")
        if (identify == "mouse") command <- suppressMarkdown(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qqPlot", reset = "QQPlot", apply="QQPlot")
    distFrame <- tkframe(optFrame)
    distVariable <- tclVar(dialog.values$initial.dist)
    normalButton <- ttkradiobutton(distFrame, variable = distVariable,
        value = "norm")
    tButton <- ttkradiobutton(distFrame, variable = distVariable,
        value = "t")
    chisqButton <- ttkradiobutton(distFrame, variable = distVariable,
        value = "chisq")
    FButton <- ttkradiobutton(distFrame, variable = distVariable,
        value = "f")
    otherButton <- ttkradiobutton(distFrame, variable = distVariable,
        value = "other")
    tDfFrame <- tkframe(distFrame)
    tDfVariable <- tclVar(dialog.values$initial.df)
    tDfField <- ttkentry(tDfFrame, width = "6", textvariable = tDfVariable)
    chisqDfFrame <- tkframe(distFrame)
    chisqDfVariable <- tclVar(dialog.values$initial.chisqdf)
    chisqDfField <- ttkentry(chisqDfFrame, width = "6", textvariable = chisqDfVariable)
    FDfFrame <- tkframe(distFrame)
    FDf1Variable <- tclVar(dialog.values$initial.fdf1)
    FDf1Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf1Variable)
    FDf2Variable <- tclVar(dialog.values$initial.fdf2)
    FDf2Field <- ttkentry(FDfFrame, width = "6", textvariable = FDf2Variable)
    otherParamsFrame <- tkframe(distFrame)
    otherParamsVariable <- tclVar(dialog.values$initial.otherparam)
    otherParamsField <- ttkentry(otherParamsFrame, width = "30",
        textvariable = otherParamsVariable)
    otherNameVariable <- tclVar(dialog.values$initial.othername)
    otherNameField <- ttkentry(otherParamsFrame, width = "10",
        textvariable = otherNameVariable)
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(labelRcmdr(distFrame, text = gettextRcmdr("Distribution"),
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), columnspan = 6, sticky = "w")
    tkgrid(normalButton, labelRcmdr(distFrame, text = gettextRcmdr("Normal")),
        sticky = "w")
    tkgrid(labelRcmdr(tDfFrame, text = gettextRcmdr("df = ")),
        tDfField, sticky = "w")
    tkgrid(tButton, labelRcmdr(distFrame, text = "t"), tDfFrame,
        sticky = "w")
    tkgrid(labelRcmdr(chisqDfFrame, text = gettextRcmdr("df = ")),
        chisqDfField, sticky = "w")
    tkgrid(chisqButton, labelRcmdr(distFrame, text = gettextRcmdr("Chi-square")),
        chisqDfFrame, sticky = "w")
    tkgrid(labelRcmdr(FDfFrame, text = gettextRcmdr("Numerator df = ")),
        FDf1Field, labelRcmdr(FDfFrame, text = gettextRcmdr("Denominator df = ")),
        FDf2Field, sticky = "w")
    tkgrid(FButton, labelRcmdr(distFrame, text = "F"), FDfFrame,
        sticky = "w")
    tkgrid(labelRcmdr(otherParamsFrame, text = gettextRcmdr("Specify: ")),
        otherNameField, labelRcmdr(otherParamsFrame, text = gettextRcmdr("Parameters: ")),
        otherParamsField, sticky = "w")
    tkgrid(otherButton, labelRcmdr(distFrame, text = gettextRcmdr("Other")),
        otherParamsFrame, sticky = "w")
    tkgrid(distFrame, sticky = "w")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(identifyPointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(identifyPointsFrame, sticky="w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

PlotMeans <- function () {
    defaults <- list(initial.groups = NULL, initial.response = NULL, initial.error.bars = "se",
        initial.level = "0.95", initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"),
        initial.main=gettextRcmdr("<auto>"), initial.tab=0)
    dialog.values <- getDialog("PlotMeans", defaults)
    initializeDialog(title = gettextRcmdr("Plot Means"), use.tabs=TRUE)
    groupBox <- variableListBox(dataTab, Factors(), title = gettextRcmdr("Factors (pick one or two)"),
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.groups, "factor"))
    responseBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.response, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Error Bars"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
        closeDialog()
        if (0 == length(groups)) {
            errorCondition(recall = PlotMeans, message = gettextRcmdr("No factors selected."))
            return()
        }
        if (2 < length(groups)) {
            errorCondition(recall = PlotMeans, message = gettextRcmdr("More than two factors selected."))
            return()
        }
        if (0 == length(response)) {
            errorCondition(recall = PlotMeans, message = gettextRcmdr("No response variable selected."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        error.bars <- tclvalue(errorBarsVariable)
        level <- if (error.bars == "conf.int")
            paste(", level=", tclvalue(levelVariable), sep = "")
        else ""
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("PlotMeans", list(initial.groups = groups, initial.response = response,
            initial.error.bars = error.bars,
            initial.level = tclvalue(levelVariable),
            initial.xlab=tclvalue(xlabVar), initial.ylab=tclvalue(ylabVar),
            initial.main=tclvalue(mainVar), initial.tab=tab))
        if (length(groups) == 1)
            doItAndPrint(paste("with(", .activeDataSet, ", plotMeans(", 
                response, ", ", groups[1],
                ", error.bars=\"", error.bars, "\"", level, xlab, ylab, main, "))",
                sep = ""))
        else {
            if (eval(parse(text = paste("with(", .activeDataSet, ", length(levels(", 
                groups[1], ")) < length(levels(", 
                groups[2], ")))", sep = ""))))
                groups <- rev(groups)
            doItAndPrint(paste("with(", .activeDataSet, ", plotMeans(", 
                response, ", ", groups[1],
                ", ", groups[2], ", error.bars=\"",
                error.bars, "\"", level, xlab, ylab, main, "))", sep = ""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    errorBarsVariable <- tclVar(dialog.values$initial.error.bars)
    seButton <- ttkradiobutton(optFrame, variable = errorBarsVariable,
        value = "se")
    sdButton <- ttkradiobutton(optFrame, variable = errorBarsVariable,
        value = "sd")
    confIntButton <- ttkradiobutton(optFrame, variable = errorBarsVariable,
        value = "conf.int")
    noneButton <- ttkradiobutton(optFrame, variable = errorBarsVariable,
        value = "none")
    levelVariable <- tclVar(dialog.values$initial.level)
    levelEntry <- ttkentry(optFrame, width = "6", textvariable = levelVariable)
    OKCancelHelp(helpSubject = "plotMeans", reset = "PlotMeans", apply = "PlotMeans")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
    #     tkgrid(labelRcmdr(optFrame, text = gettextRcmdr("Error Bars"),
    #                       fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    tkgrid(seButton, labelRcmdr(optFrame, text = gettextRcmdr("Standard errors")),
        sticky = "w")
    tkgrid(sdButton, labelRcmdr(optFrame, text = gettextRcmdr("Standard deviations")),
        sticky = "w")
    tkgrid(confIntButton, labelRcmdr(optFrame, text = gettextRcmdr("Confidence intervals")),
        labelRcmdr(optFrame, text = gettextRcmdr("   Level of confidence:")),
        levelEntry, sticky = "w")
    tkgrid(noneButton, labelRcmdr(optFrame, text = gettextRcmdr("No error bars")),
        sticky = "w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

Scatter3D <- function () {
    use.rgl <- getOption("Rcmdr")$use.rgl
    if (is.null(use.rgl) || use.rgl) {
        Library("rgl")
        Library("mgcv")
    }
    putRcmdr("rgl.command", TRUE)
    defaults <- list (initial.x = NULL, initial.y = NULL, initial.scales = 1, initial.grid = 1, 
        initial.resids = 0, initial.lin = 0, initial.quad = 0, initial.nonpar = 0, 
        initial.additive = 0, initial.ellips = 0, initial.dfNonpar = gettextRcmdr("<auto>"), 
        initial.dfAdd = gettextRcmdr("<auto>"), initial.bg = "white",
        initialGroup=NULL, initial.lines.by.group=0, initial.identify="not", initial.id.n="2",
        initial.tab=0)
    dialog.values <- getDialog ("Scatter3D", defaults)
    initial.group <- dialog.values$initial.group
    .linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    initializeDialog(title = gettextRcmdr("3D Scatterplot"), use.tabs=TRUE)
    variablesFrame <- tkframe(dataTab)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("Explanatory variables (pick two)"), 
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    yBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("Response variable (pick one)"), 
        initialSelection = varPosn (dialog.values$initial.y, "numeric"))
    surfacesFrame <- tkframe(optionsTab)
    axisScales <- tclVar(dialog.values$initial.scales)
    axisScalesCheckBox <- ttkcheckbutton(surfacesFrame, variable = axisScales)
    gridLines <- tclVar(dialog.values$initial.grid)
    gridLinesCheckBox <- ttkcheckbutton(surfacesFrame, variable = gridLines)
    squaredResiduals <- tclVar(dialog.values$initial.resids)
    squaredResidualsCheckBox <- ttkcheckbutton(surfacesFrame, 
        variable = squaredResiduals)
    linearLSSurface <- tclVar(dialog.values$initial.lin)
    linearLSCheckBox <- ttkcheckbutton(surfacesFrame, variable = linearLSSurface)
    quadLSSurface <- tclVar(dialog.values$initial.quad)
    quadLSCheckBox <- ttkcheckbutton(surfacesFrame, variable = quadLSSurface)
    nonparSurface <- tclVar(dialog.values$initial.nonpar)
    nonparCheckBox <- ttkcheckbutton(surfacesFrame, variable = nonparSurface)
    dfNonparVariable <- tclVar(dialog.values$initial.dfNonpar)
    dfNonparField <- ttkentry(surfacesFrame, width = "6", textvariable = dfNonparVariable)
    additiveSurface <- tclVar(dialog.values$initial.additive)
    additiveCheckBox <- ttkcheckbutton(surfacesFrame, variable = additiveSurface)
    dfAddVariable <- tclVar(dialog.values$initial.dfAdd)
    dfAddField <- ttkentry(surfacesFrame, width = "6", textvariable = dfAddVariable)
    ellipsoid <- tclVar(dialog.values$initial.ellips)
    ellipsoidCheckBox <- ttkcheckbutton(surfacesFrame, variable = ellipsoid)
    bgFrame <- tkframe(optionsTab)
    bgVariable <- tclVar(dialog.values$initial.bg)
    whiteButton <- ttkradiobutton(bgFrame, variable = bgVariable, 
        value = "white")
    blackButton <- ttkradiobutton(bgFrame, variable = bgVariable, 
        value = "black")
    idFrame <- tkframe(optionsTab)
    radioButtons(window=idFrame, name = "identify", buttons = c("auto", "mouse", 
        "not"), labels = gettextRcmdr(c("Automatically", 
            "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"), 
        initialValue = dialog.values$initial.identify)
    id.n.Var <- tclVar(dialog.values$initial.id.n) 
    npointsSpinner <- tkspinbox(idFrame, from=1, to=10, width=2, textvariable=id.n.Var)    
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        scales <- tclvalue(axisScales)
        grid <- tclvalue(gridLines)
        resids <- tclvalue(squaredResiduals)
        lin <- tclvalue(linearLSSurface)
        quad <- tclvalue(quadLSSurface)
        nonpar <- tclvalue(nonparSurface)
        additive <- tclvalue(additiveSurface)
        ellips <- tclvalue(ellipsoid) 
        dfNonpar <- tclvalue(dfNonparVariable)
        dfAdd <- tclvalue(dfAddVariable)
        bg <- tclvalue(bgVariable)
        identify <- tclvalue(identifyVariable)
        id.n <- tclvalue(id.n.Var)
        identify.text <- switch(identify,
            auto = paste(", id.method='mahal', id.n =", id.n),
            mouse = ", id.method='identify'",
            not = "")
        closeDialog()
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = scatterPlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        if (length(y) == 0) {
            errorCondition(recall = Scatter3D, message = gettextRcmdr("You must select a response variable."))
            return()
        }
        if (2 != length(x)) {
            errorCondition(recall = Scatter3D, message = gettextRcmdr("You must select 2 explanatory variables."))
            return()
        }
        if (is.element(y, x)) {
            errorCondition(recall = Scatter3D, message = gettextRcmdr("Response and explanatory variables must be different."))
            return()
        }
        putDialog ("Scatter3D", list(initial.x = x, initial.y = y, initial.scales = scales, initial.grid = grid, 
            initial.resids = resids, initial.lin = lin, initial.quad = quad, initial.nonpar = nonpar, 
            initial.additive = additive, initial.ellips = ellips, initial.dfNonpar = dfNonpar, 
            initial.dfAdd = dfAdd, initial.bg = bg, 
            initial.group=if (.groups == FALSE) NULL else .groups,
            initial.lines.by.group=if (.linesByGroup) 1 else 0,
            initial.identify=identify, initial.id.n=id.n, initial.tab=tab))
        scales <- if (tclvalue(axisScales) == 1) 
            "TRUE"
        else "FALSE"
        grid <- if (tclvalue(gridLines) == 1) 
            "TRUE"
        else "FALSE"
        resids <- if (tclvalue(squaredResiduals) == 1) 
            ", residuals=\"squares\""
        else ", residuals=TRUE"
        lin <- if (tclvalue(linearLSSurface) == 1) 
            "\"linear\""
        quad <- if (tclvalue(quadLSSurface) == 1) 
            "\"quadratic\""
        nonpar <- if (tclvalue(nonparSurface) == 1) 
            "\"smooth\""
        additive <- if (tclvalue(additiveSurface) == 1) 
            "\"additive\""
        surfaces <- c(lin, quad, nonpar, additive)
        nsurfaces <- length(surfaces)
        if (nsurfaces > 1) 
            resids <- ""
        ellips <- if (tclvalue(ellipsoid) == 1) 
            "TRUE"
        else "FALSE"
        opts <- options(warn = -1)
        dfNonpar <- if (dfNonpar == gettextRcmdr("<auto>")) 
            ""
        else paste(", df.smooth=", as.numeric(dfNonpar), sep = "")
        dfAdd <- if (dfAdd == gettextRcmdr("<auto>")) 
            ""
        else paste(", df.additive=", as.numeric(dfAdd), sep = "")
        options(opts)
        fit <- if (nsurfaces == 0) 
            ", surface=FALSE"
        else if (nsurfaces == 1) 
            paste(", fit=", surfaces, sep = "")
        else paste(", fit=c(", paste(surfaces, collapse = ","), 
            ")", sep = "")
        .activeDataSet <- ActiveDataSet()
        if (.groups != FALSE) {
            groups <- paste(", groups=", .activeDataSet, "$", 
                .groups, sep = "")
            parallel <- paste(", parallel=", .linesByGroup, sep = "")
        }
        else parallel <- groups <- ""
        if (identify == "mouse"){
            RcmdrTkmessageBox(title="Identify Points",
                message=gettextRcmdr("Drag right mouse button to identify points,\nclick right button to exit."),
                icon="info", type="ok")
        }
        if (.groups == FALSE) {
            command <- paste("scatter3d(", y, "~", x[1], "+", x[2], ", data=", .activeDataSet, 
                fit, resids, dfNonpar, dfAdd, 
                parallel, ", bg=\"", bg, "\", axis.scales=", scales, 
                ", grid=", grid, ", ellipsoid=", ellips, identify.text,
                ")", sep = "")
            if (identify == "mouse") command <- suppressMarkdown(command)
            doItAndPrint(command)
        }
        else {
            command <- paste("scatter3d(", y, "~", x[1], "+", x[2], "|", .groups, ", data=", .activeDataSet, 
                fit, resids, dfNonpar, dfAdd, 
                parallel, ", bg=\"", bg, "\", axis.scales=", scales, 
                ", grid=", grid, ", ellipsoid=", ellips, identify.text,
                ")", sep = "")
            if (identify == "mouse") command <- suppressMarkdown(command)
            doItAndPrint(command)
        }
        putRcmdr("rgl", TRUE)
        command <- paste("Identify3d(", .activeDataSet, "$", 
            x[1], ", ", .activeDataSet, "$", y, ", ", .activeDataSet, 
            "$", x[2], groups, ", axis.scales=", scales, ", labels=row.names(", 
            .activeDataSet, "))", sep = "")
        putRcmdr("Identify3d", command)
        .Tcl("update")
        activateMenus()
        tkfocus(CommanderWindow())
        rgl.bringtotop()
    }
    groupsBox(Scatter3D, plotLinesByGroup = TRUE, plotLinesByGroupsText = gettextRcmdr("Parallel regression surfaces"),
        initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group),
        window=dataTab)
    OKCancelHelp(helpSubject = "Scatter3DDialog", reset = "Scatter3D", apply = "Scatter3D")
    tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text = "  "), 
        getFrame(yBox), sticky = "nw")
    tkgrid(variablesFrame, sticky = "nw")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Show axis scales")), 
        axisScalesCheckBox, sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Show surface grid lines")), 
        gridLinesCheckBox, sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Show squared residuals")), 
        squaredResidualsCheckBox, sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Surfaces to Fit"), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Linear least-squares")), 
        linearLSCheckBox, sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Quadratic least-squares")), 
        quadLSCheckBox, sticky = "w")
    dfLabel <- labelRcmdr(surfacesFrame, text = gettextRcmdr("df = "))
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Smooth regression")), 
        nonparCheckBox, dfLabel, dfNonparField, sticky = "w")
    tkgrid.configure(dfLabel, sticky = "e")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Additive regression")), 
        additiveCheckBox, labelRcmdr(surfacesFrame, text = gettextRcmdr("df(each term) = ")), 
        dfAddField, sticky = "w")
    tkgrid(labelRcmdr(surfacesFrame, text = gettextRcmdr("Plot 50% concentration ellipsoid")), 
        ellipsoidCheckBox, sticky = "w")
    tkgrid(surfacesFrame, sticky = "w")
    tkgrid(labelRcmdr(bgFrame, text = gettextRcmdr("Background Color"), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w", columnspan = 2)
    tkgrid(labelRcmdr(bgFrame, text = gettextRcmdr("Black")), 
        blackButton, sticky = "w")
    tkgrid(labelRcmdr(bgFrame, text = gettextRcmdr("White")), 
        whiteButton, sticky = "w")
    tkgrid(bgFrame, sticky = "w")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(idFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(idFrame, sticky="w")
    tkgrid(groupsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

Identify3D <- function(){
    if (0 == rgl.cur()) {
        Message(message=gettextRcmdr("There is no current RGL graphics device."),
            type="error")
        return()
    }
    RcmdrTkmessageBox(title="Identify Points",
        message=gettextRcmdr("Drag right mouse button to identify points,\nclick right button to exit."),
        icon="info", type="ok")
    command <- getRcmdr("Identify3d")
    command <- suppressMarkdown(command)
    doItAndPrint(command)
}

saveBitmap <- function () {
    env <- environment()
    updateWidth <- function(...){
        if (tclvalue(aspectVariable) == "1"){
            tclvalue(heightVariable) <- round(aspect*as.numeric(tclvalue(widthVariable)))
        }
    }
    updateHeight <- function(...){
        if (tclvalue(aspectVariable) == "1"){
            tclvalue(widthVariable) <- round((1/aspect)*as.numeric(tclvalue(heightVariable)))
        }
    }
    updateSize <- function(...){
        units <- tclvalue(unitsVariable)
        size <- dev.size(units=units)
        if (units == "in") {
            wmin <- min(3, size[1])
            wmax <- max(10, size[1])
            hmin <- min(3, size[2])
            hmax <- max(10, size[2])
            rmin <- 50
            rmax <- 300
            res <- if (tclvalue(resVariable) == "72") 72 else round(2.54*as.numeric(tclvalue(resVariable)))
        }
        else if (units == "cm") {
            wmin <- min(8, size[1])
            wmax <- max(25, size[1])
            hmin <- min(8, size[2])
            hmax <- max(25, size[2])
            rmin <- 20
            rmax <- 120
            res <- round(as.numeric(tclvalue(resVariable))/2.54)
        }
        else {
            wmin <- min(200, size[1])
            wmax <- max(1000, size[1])
            hmin <- min(200, size[2])
            hmax <- max(1000, size[2])
            rmin <- 50
            rmax <- 300
            res <- 72
        }
        tkconfigure(widthSlider, from = wmin, to = wmax)
        tkconfigure(heightSlider,  from = hmin, to = hmax)
        tkconfigure(wlabel, text = paste(gettextRcmdr(c("Width", " (", all.units[units], ")")), collapse=""))
        tkconfigure(hlabel, text = paste(gettextRcmdr(c("Height",  " (", all.units[units], ")")), collapse=""))
        tkconfigure(rlabel, text = paste(gettextRcmdr(c("Resolution (pixels/", unit[units], ")")), collapse=""))
        tkconfigure(resSlider, from=rmin, to=rmax, state = if (tclvalue(unitsVariable) == "px") "disabled" else "normal")
        tkconfigure(disabled, text = if (units == "px") gettextRcmdr("[disabled]") else "")
        tclvalue(widthVariable) <- size[1]
        tclvalue(heightVariable) <- size[2]
        tclvalue(resVariable) <- res
    }
    all.units <- c("inches", "cm", "pixels")
    names(all.units) <- c("in", "cm", "px")
    unit <- c("inch", "cm", "inch")
    names(unit) <- c("in", "cm", "px")
    if (1 == dev.cur()) {
        Message(gettextRcmdr("There is no current graphics device to save."), 
            type = "error")
        return()
    }
    defaults <- list (initial.type = "png", initial.pointsize=12, initial.units="in", initial.res = 72)
    dialog.values <- getDialog ("saveBitmap", defaults)
    units <- dialog.values$initial.units
    size <- dev.size(units=units)
    aspect <- size[2]/size[1]
    if (units == "in") {
        wmin <- min(3, size[1])
        wmax <- max(10, size[1])
        hmin <- min(3, size[2])
        hmax <- max(10, size[2])
        rmin <- 50
        rmax <- 300
        res <- dialog.values$initial.res
    }
    else if (units == "cm") {
        wmin <- min(8, size[1])
        wmax <- max(25, size[1])
        hmin <- min(8, size[2])
        hmax <- max(25, size[2])
        rmin <- 20
        rmax <- 120
        res <- dialog.values$initial.res
    }
    else {
        wmin <- min(200, size[1])
        wmax <- max(1000, size[1])
        hmin <- min(200, size[2])
        hmax <- max(1000, size[2])
        rmin <- 50
        rmax <- 300
        res <- 72
    }
    initializeDialog(title = gettextRcmdr("Save Graph as Bitmap"))
    radioButtons(name = "filetype", buttons = c("png", "jpeg"), 
        labels = c("PNG", "JPEG"), title = gettextRcmdr("Graphics File Type"),
        initialValue = dialog.values$initial.type)
    radioButtons(name = "units", buttons = c("in", "cm", "px"), 
        labels = gettextRcmdr(c("inches", "cm", "pixels")), title = gettextRcmdr("Units"),
        initialValue = dialog.values$initial.units, command=updateSize)
    sliderFrame <- tkframe(top)
    widthVariable <- tclVar(size[1])
    widthSlider <- tkscale(sliderFrame, from = wmin, to = wmax, 
        showvalue = TRUE, variable = widthVariable, resolution = 1, 
        orient = "horizontal", command=updateWidth)
    heightVariable <- tclVar(size[2])
    heightSlider <- tkscale(sliderFrame, from = hmin, to = hmax, 
        showvalue = TRUE, variable = heightVariable, resolution = 1, 
        orient = "horizontal", command=updateHeight)
    pointSizeVariable <- tclVar(dialog.values$initial.pointsize)
    pointSizeSlider <- tkscale(sliderFrame, from = 6, to = 16, 
        showvalue = TRUE, variable = pointSizeVariable, resolution = 1, 
        orient = "horizontal")
    resVariable <- tclVar(res)
    resSlider <- tkscale(sliderFrame, from = rmin, to = rmax, 
        showvalue = TRUE, variable = resVariable, resolution = 1, 
        orient = "horizontal")
    tkconfigure(resSlider,  state = if (tclvalue(unitsVariable) == "px") "disabled" else "normal")
    aspectVariable <- tclVar("1")
    aspectFrame <- tkframe(top)
    aspectCheckBox <- ttkcheckbutton(aspectFrame, variable = aspectVariable)
    onOK <- function() {
        closeDialog()
        width <- tclvalue(widthVariable)
        height <- tclvalue(heightVariable)
        type <- tclvalue(filetypeVariable)
        pointsize <- tclvalue(pointSizeVariable)
        units <- tclvalue(unitsVariable)
        res <- tclvalue(resVariable)
        putDialog ("saveBitmap", list (initial.type = type, initial.pointsize = pointsize, initial.units=units, initial.res=res))
        if (type == "png") {
            ext <- "png"
            filetypes <- gettextRcmdr("{\"All Files\" {\"*\"}} {\"PNG Files\" {\".png\" \".PNG\"}}")
            initial <- "RGraph.png"
        }
        else {
            ext <- "jpg"
            filetypes <- gettextRcmdr("{\"All Files\" {\"*\"}} {\"JPEG Files\" {\".jpg\" \".JPG\" \".jpeg\" \".JPEG\"}}")
            initial <- "RGraph.jpg"
        }
        filename <- tclvalue(tkgetSaveFile(filetypes = filetypes, 
            defaultextension = ext, initialfile = initial, parent = CommanderWindow()))
        if (filename == "") 
            return()
        command <- paste("dev.print(", type, ", filename=\"", 
            filename, "\", width=", width, ", height=", height, ", pointsize=", pointsize, ', units="', units, 
            if(units == "px") '")' else paste('", res=', res, ')', sep=""), sep = "")
        doItAndPrint(command, rmd=FALSE)
        Message(paste(gettextRcmdr("Graph saved to file"), filename), 
            type = "note")
    }
    OKCancelHelp(helpSubject = "png", reset = "saveBitmap")
    tkgrid(filetypeFrame, sticky = "w")
    tkgrid(unitsFrame, stick="w")
    tkgrid(labelRcmdr(aspectFrame, text = gettextRcmdr("Fixed aspect ratio (height:width)")),
        aspectCheckBox, sticky="w")
    tkgrid(aspectFrame, sticky="w")
    tkgrid(wlabel <- labelRcmdr(sliderFrame, text = paste(gettextRcmdr(c("Width", " (", all.units[units], ")")), collapse="")), 
        widthSlider, sticky = "sw")
    tkgrid(hlabel <- labelRcmdr(sliderFrame, text = paste(gettextRcmdr(c("Height",  " (", all.units[units], ")")), collapse="")), 
        heightSlider, sticky = "sw")
    tkgrid(rlabel <- labelRcmdr(sliderFrame, text = paste(gettextRcmdr(c("Resolution", "(", "pixels", "/", unit[units], ")")), collapse="")), 
        resSlider, 
        disabled <- labelRcmdr(sliderFrame, text = if (units == "px") gettextRcmdr("[disabled]") else ""),
        sticky = "sw")
    tkgrid(labelRcmdr(sliderFrame, text = gettextRcmdr("Text size (points)")), 
        pointSizeSlider, sticky = "sw")
    tkgrid(sliderFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

savePDF <- function () {
    updateWidth <- function(...){
        if (tclvalue(aspectVariable) == "1"){
            tclvalue(heightVariable) <- round(aspect*as.numeric(tclvalue(widthVariable)), 1)
        }
    }
    updateHeight <- function(...){
        if (tclvalue(aspectVariable) == "1"){
            tclvalue(widthVariable) <- round((1/aspect)*as.numeric(tclvalue(heightVariable)), 1)
        }
    }
    updateSize <- function(...){
        units <- tclvalue(unitsVariable)
        size <- dev.size(units=units)
        if (units == "in") {
            wmin <- min(3, size[1])
            wmax <- max(10, size[1])
            hmin <- min(3, size[2])
            hmax <- max(10, size[2])
        }
        else {
            wmin <- min(8, size[1])
            wmax <- max(25, size[1])
            hmin <- min(8, size[2])
            hmax <- max(25, size[2])
        }
        tkconfigure(widthSlider, from = wmin, to = wmax)
        tkconfigure(heightSlider,  from = hmin, to = hmax)
        tkconfigure(wlabel, text = paste(gettextRcmdr(c("Width", " (", all.units[units], ")")), collapse=""))
        tkconfigure(hlabel, text = paste(gettextRcmdr(c("Height",  " (", all.units[units], ")")), collapse=""))
        tclvalue(widthVariable) <- size[1]
        tclvalue(heightVariable) <- size[2]
    }
    all.units <- c("inches", "cm")
    names(all.units) <- c("in", "cm")
    if (1 == dev.cur()) {
        Message(gettextRcmdr("There is no current graphics device to save."), 
            type = "error")
        return()
    }
    defaults <- list (initial.type = "pdf", initial.pointsize = 12, initial.units="in")
    dialog.values <- getDialog ("savePDF", defaults)
    units <- dialog.values$initial.units
    size <- dev.size(units=units)
    aspect <- size[2]/size[1]
    size <- round(size, 1)
    if (units == "in") {
        wmin <- min(3, size[1])
        wmax <- max(10, size[1])
        hmin <- min(3, size[2])
        hmax <- max(10, size[2])
    }
    else {
        wmin <- min(8, size[1])
        wmax <- max(25, size[1])
        hmin <- min(8, size[2])
        hmax <- max(25, size[2])
    }
    initializeDialog(title = gettextRcmdr("Save Graph as PDF/Postscript"))
    radioButtons(name = "filetype", buttons = c("pdf", "postscript", 
        "eps"), labels = gettextRcmdr(c("PDF", "Postscript", 
            "Encapsulated Postscript")), title = gettextRcmdr("Graphics File Type"), 
        initialValue = dialog.values$initial.type)
    radioButtons(name = "units", buttons = c("in", "cm"), 
        labels = gettextRcmdr(c("inches", "cm")), title = gettextRcmdr("Units"),
        initialValue = dialog.values$initial.units, command=updateSize)
    aspectVariable <- tclVar("1")
    aspectFrame <- tkframe(top)
    aspectCheckBox <- ttkcheckbutton(aspectFrame, variable = aspectVariable)
    sliderFrame <- tkframe(top)
    widthVariable <- tclVar(size[1])
    widthSlider <- tkscale(sliderFrame, from = wmin, to = wmax, 
        showvalue = TRUE, 
        variable = widthVariable, resolution = 0.1, orient = "horizontal", 
        command=updateWidth)
    heightVariable <- tclVar(size[2])
    heightSlider <- tkscale(sliderFrame, from = hmin, to = hmax, 
        showvalue = TRUE, 
        variable = heightVariable, resolution = 0.1, orient = "horizontal",
        command=updateHeight)
    pointSizeVariable <- tclVar(dialog.values$initial.pointsize)
    pointSizeSlider <- tkscale(sliderFrame, from = 6, to = 16, 
        showvalue = TRUE, variable = pointSizeVariable, resolution = 1, 
        orient = "horizontal")
    onOK <- function() {
        closeDialog()
        width <- tclvalue(widthVariable)
        height <- tclvalue(heightVariable)
        type <- tclvalue(filetypeVariable)
        units <- tclvalue(unitsVariable)
        pointsize <- tclvalue(pointSizeVariable)
        putDialog ("savePDF", list (initial.type = type, initial.pointsize = pointsize, initial.units=units))
        if (units == "cm") {
            width <- round(as.numeric(width)/2.54, 1)
            height <- round(as.numeric(height)/2.54, 1)
        } 
        if (type == "pdf") {
            ext <- "pdf"
            filetypes <- gettextRcmdr("{\"All Files\" {\"*\"}} {\"PDF Files\" {\".pdf\" \".PDF\"}}")
            initial <- "RGraph.pdf"
        }
        else if (type == "postscript") {
            ext <- "ps"
            filetypes <- gettextRcmdr("{\"All Files\" {\"*\"}} {\"Postscript Files\" {\".ps\" \".PS\"}}")
            initial <- "RGraph.ps"
        }
        else {
            ext <- "eps"
            filetypes <- gettextRcmdr("{\"All Files\" {\"*\"}} {\"Encapsulated Postscript Files\" {\".eps\" \".EPS\"}}")
            initial <- "RGraph.eps"
        }
        filename <- tclvalue(tkgetSaveFile(filetypes = filetypes, 
            defaultextension = ext, initialfile = initial, parent = CommanderWindow()))
        if (filename == "") 
            return()
        command <- if (type == "eps") 
            paste("dev.copy2eps(file=\"", filename, "\", width=", 
                width, ", height=", height, ", pointsize=", pointsize, 
                ")", sep = "")
        else paste("dev.print(", type, ", file=\"", filename, 
            "\", width=", width, ", height=", height, ", pointsize=", 
            pointsize, ")", sep = "")
        doItAndPrint(command, rmd=FALSE)
        Message(paste(gettextRcmdr("Graph saved to file"), filename), 
            type = "note")
    }
    OKCancelHelp(helpSubject = "pdf", reset = "savePDF")
    tkgrid(filetypeFrame, sticky = "w")
    tkgrid(unitsFrame, stick="w")
    tkgrid(labelRcmdr(aspectFrame, text = gettextRcmdr("Fixed aspect ratio (height:width)")),
        aspectCheckBox, sticky="w")
    tkgrid(aspectFrame, sticky="w")
    tkgrid(wlabel <- labelRcmdr(sliderFrame, text = paste(gettextRcmdr(c("Width", " (", all.units[units], ")")), collapse="")), 
        widthSlider, sticky = "sw")
    tkgrid(hlabel <- labelRcmdr(sliderFrame, text = paste(gettextRcmdr(c("Height",  " (", all.units[units], ")")), collapse="")), 
        heightSlider, sticky = "sw")
    tkgrid(labelRcmdr(sliderFrame, text = gettextRcmdr("Text size (points)")), 
        pointSizeSlider, sticky = "sw")
    tkgrid(sliderFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

saveRglGraph <- function(){
    if (0 == rgl.cur()) {
        Message(message=gettextRcmdr("There is no current RGL graphics device to save."),
            type="error")
        return()
    }
    ext <- "png"
    filetypes <- gettextRcmdr('{"All Files" {"*"}} {"PNG Files" {".png" ".PNG"}}')
    initial <- "RGLGraph.png"
    filename <- tclvalue(tkgetSaveFile(filetypes=filetypes,
        defaultextension=ext,
        initialfile=initial,
        parent=CommanderWindow()))
    if (filename == "") return()
    command <- paste('rgl.snapshot("', filename, '")', sep="")
    doItAndPrint(command, rmd=FALSE)
    Message(paste(gettextRcmdr("Graph saved to file"), filename), type="note")
}

## The following function by Richard Heiberger, with small modifications by J. Fox
## with more modifications by Richard Heiberger.
## 2008-01-03 added conditions, layout, and multiple colors
## 2012-08-19 rmh added memory to the dialogs, using John Fox's getDialog and putDialog functions
## 2013-06-19 J. Fox added Data and Options tabs, Apply button

Xyplot <- function() {
    Library("lattice")
    defaults <- list(initial.predictor = NULL, initial.response = NULL,
        initial.auto.key = 1, initial.outer = 0,
        initial.x.relation = "same", initial.y.relation = "same",
        initial.layoutColumns = "", initial.layoutRows = "",
        initial.conditions = FALSE,
        initial.groups = FALSE,
        initial.points = 1, initial.lines = 0, initial.tab=0,
        initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"),
        initial.main=gettextRcmdr("<auto>"))
    dialog.values <- getDialog("Xyplot", defaults)
    initializeDialog(title=gettextRcmdr("XY Conditioning Plot"), use.tabs=TRUE)
    predictorFrame <- tkframe(dataTab)
    predictorBox <-
        variableListBox(predictorFrame, Numeric(),
            title=gettextRcmdr("Explanatory variables (pick one or more)"),
            selectmode="multiple",
            initialSelection = varPosn (dialog.values$initial.predictor, "numeric"))
    responseBox <- variableListBox(predictorFrame, Numeric(),
        title=gettextRcmdr("Response variables (pick one or more)"),
        selectmode="multiple",
        initialSelection = varPosn (dialog.values$initial.response, "numeric"))
    cgFrame <- tkframe(dataTab)
    conditions.if <-
        length(dialog.values$initial.conditions) == 1 &&
        dialog.values$initial.conditions == FALSE
    conditionsBox <- variableListBox(cgFrame, Factors(),
        title=gettextRcmdr("Conditions '|' (pick zero or more)"),
        selectmode="multiple",
        initialSelection=if (conditions.if) FALSE else
            varPosn(dialog.values$initial.conditions, "factor"))
    groups.if <-
        length(dialog.values$initial.groups) == 1 &&
        dialog.values$initial.groups == FALSE
    groupsBox <- variableListBox(cgFrame, Factors(),
        title=gettextRcmdr("Groups 'groups=' (pick zero or more)"),
        selectmode="multiple",
        initialSelection=if (groups.if) FALSE else
            varPosn(dialog.values$initial.groups, "factor"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Options"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    checkBoxes(window = optFrame, frame="otherFrame",
        boxes=c("auto.key", "outer"),
        initialValues=c(dialog.values$initial.auto.key, dialog.values$initial.outer),
        labels=gettextRcmdr(c("Automatically draw key",
            "Different panels for different y ~ x combinations")))
    relationFrame <- tkframe(optFrame)
    radioButtons(window=relationFrame,
        name="x.relation",
        buttons=c("same", "free", "sliced"),
        labels=gettextRcmdr(c("Identical", "Free", "Same range")),
        title=gettextRcmdr("X-Axis Scales in Different Panels"),
        initialValue = dialog.values$initial.x.relation)
    radioButtons(window=relationFrame,
        name="y.relation",
        buttons=c("same", "free", "sliced"),
        labels=gettextRcmdr(c("Identical", "Free", "Same range")),
        title=gettextRcmdr("Y-Axis Scales in Different Panels"),
        initialValue = dialog.values$initial.y.relation)
    
    scalarsFrame <- tkframe(optFrame)
    layoutColumnsVar <- tclVar(dialog.values$initial.layoutColumns)
    layoutColumnsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutColumnsVar)
    layoutRowsVar <- tclVar(dialog.values$initial.layoutRows)
    layoutRowsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutRowsVar)
    
    checkBoxes(window = optFrame, frame="typeFrame",
        boxes=c("points", "lines"),
        initialValues=c(dialog.values$initial.points, dialog.values$initial.lines),
        labels=gettextRcmdr(c("Points", "Lines")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        predictor <- getSelection(predictorBox)
        response <- getSelection(responseBox)
        conditions <- getSelection(conditionsBox)
        groups <- getSelection(groupsBox)
        closeDialog()
        
        if (0 == length(response)) {
            errorCondition(recall=Xyplot,
                message=gettextRcmdr("At least one response variable must be selected."))
            return()
        }
        if (0 == length(predictor)) {
            errorCondition(recall=Xyplot,
                message=gettextRcmdr("At least one explanatory variable must be selected."))
            return()
        }
        auto.key <- ("1" == tclvalue(auto.keyVariable))
        outer    <- ("1" == tclvalue(outerVariable))
        x.relation <- as.character(tclvalue(x.relationVariable))
        y.relation <- as.character(tclvalue(y.relationVariable))
        
        layoutColumns  <- as.numeric(tclvalue(layoutColumnsVar))
        layoutRows     <- as.numeric(tclvalue(layoutRowsVar))
        
        points <- ("1" == tclvalue(pointsVariable))
        lines  <- ("1" == tclvalue(linesVariable))
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("Xyplot", list(initial.predictor = predictor, initial.response = response,
            initial.auto.key = auto.key, initial.outer = outer,
            initial.x.relation = x.relation,
            initial.y.relation = y.relation,
            initial.layoutColumns = tclvalue(layoutColumnsVar),
            initial.layoutRows = tclvalue(layoutRowsVar),
            initial.conditions = if (length(conditions) != 0) conditions else FALSE,
            initial.groups = if (length(groups) != 0) groups else FALSE,
            initial.points = points,
            initial.lines = lines,
            initial.tab=tab, initial.xlab=tclvalue(xlabVar),
            initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar)))
        
        layout.command <- ""
        number.na <- is.na(layoutColumns) + is.na(layoutRows)
        
        if (number.na==1) {
            errorCondition(recall=Xyplot,
                message=gettextRcmdr("Both or neither layout values must be numbers."))
            return()
        }
        if (number.na==0) layout.command <- deparse(c(layoutColumns, layoutRows))
        
        .activeDataSet <- ActiveDataSet()
        
        
        
        conditions.command <-
            if (length(conditions) == 0) {
                if (outer) {
                    if (layout.command == "")
                        paste(", layout=c(",
                            length(predictor),
                            ",",
                            length(response),
                            ")")
                    else
                        paste(", layout=", layout.command, sep="")
                }
                else
                    if (layout.command != "")
                        paste(", layout=", layout.command, sep="")
            }
        else {  ## (length(conditions) > 0)
            if (outer) {
                condition.levels <- prod(sapply(conditions, d.f=get(.activeDataSet),
                    function(g, d.f) length(levels(d.f[[g]]))))
                if (layout.command != "")
                    paste(", layout=", layout.command, sep="")
                else
                    paste(", layout=c(",
                        condition.levels,
                        "*",
                        length(predictor),
                        ",",
                        length(response),
                        ")",
                        ## ", between=list(x=c(0,0, 1, 0,0), y=1)",
                        ", between=list(x=c(",
                        paste(rep(c(rep(0, condition.levels-1), 1),
                            length=condition.levels*length(predictor)-1),
                            collapse=","),
                        "), y=1)")
            }
            else
                if (layout.command != "")
                    paste(", layout=", layout.command, sep="")
        }
        
        
        groups.command <- switch(as.character(length(groups)),
            "0"="",
            "1"=paste(", groups=", groups, sep=""),
            paste(", groups=interaction(",
                paste(groups, collapse=","),
                ")", sep=""))
        
        if(!(points || lines)) {
            errorCondition(recall=Xyplot,
                message=gettextRcmdr("Choose at least one of points or lines."))
            return()
        }
        
        type.command <- paste(", type=",
            deparse(c("p"[points], "l"[lines])),
            sep="")
        
        xyplot.command <- paste("xyplot(",
            paste(response, collapse=" + "),
            " ~ ",
            paste(predictor, collapse=" + "),
            if (length(conditions) > 0)
                paste(" |",
                    paste(conditions, collapse=" + ")
                ) else "",
            if (outer) ", outer=TRUE",
            conditions.command,
            groups.command,
            type.command,
            ", pch=16",
            if (auto.key) ", auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16)" else "",
            paste(", scales=list(x=list(relation='",
                x.relation,
                "'), y=list(relation='",
                y.relation,
                "'))", sep=""),
            ", data=", .activeDataSet, xlab, ylab, main, ")", sep="")
        doItAndPrint(xyplot.command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="xyplot", reset = "Xyplot", apply = "Xyplot")
    tkgrid(getFrame(predictorBox), getFrame(responseBox),
        columnspan=1, sticky="w")
    tkgrid(predictorFrame, sticky="w")
    tkgrid(getFrame(conditionsBox),
        tklabel(cgFrame, text=gettextRcmdr("           ")),
        getFrame(groupsBox),
        columnspan=1, sticky="w")
    tkgrid(cgFrame, sticky="w")
    tkgrid(tklabel(optFrame, text=gettextRcmdr("Other Options"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(otherFrame, sticky="w")
    tkgrid(tklabel(optFrame, text=gettextRcmdr("Plot Type (one or both)"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(typeFrame, sticky="w")
    
    tkgrid(x.relationFrame, y.relationFrame, columnspan=2, sticky="w")
    tkgrid(relationFrame, sticky="w")
    
    tkgrid(tklabel(optFrame, text=gettextRcmdr("Layout"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"),
        sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of columns:")), layoutColumnsEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of rows:")), layoutRowsEntry, sticky="w")
    tkgrid(scalarsFrame, sticky="w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

# set the colour palette

setPalette <- function() {
    cval <- function(x,y) -sum((x-y)^2)
    contrasting <- function(x)
        optim(rep(127, 3),cval,lower=0,upper=255,method="L-BFGS-B",y=x)$par
    # the following local function from Thomas Lumley via r-help
    convert <- function (color){
        rgb <- col2rgb(color)/255
        L <- c(0.2, 0.6, 0) %*% rgb
        ifelse(L >= 0.2, "#000060", "#FFFFA0")
    }
    env <- environment()
    pal <- palette()
    pickColor <- function(initialcolor, parent){
        newcolor <- tclvalue(.Tcl(paste("tk_chooseColor", .Tcl.args(title = "Select a Color",
            initialcolor=initialcolor, parent=parent))))
        if (newcolor == "") initialcolor else newcolor
    }
    initializeDialog(title=gettextRcmdr("Set Color Palette"))
    hexcolor <- colorConverter(toXYZ = function(hex,...) {
        rgb <- t(col2rgb(hex))/255
        colorspaces$sRGB$toXYZ(rgb,...) },
        fromXYZ = function(xyz,...) {
            rgb <- colorspaces$sRGB$fromXYZ(xyz,..)
            rgb <- round(rgb,5)
            if (min(rgb) < 0 || max(rgb) > 1) as.character(NA)
            else rgb(rgb[1],rgb[2],rgb[3])},
        white = "D65", name = "#rrggbb")
    cols <- t(col2rgb(pal))
    hex <- convertColor(cols, from="sRGB", to=hexcolor, scale.in=255, scale.out=NULL)
    for (i in 1:8) assign(paste("hex", i, sep="."), hex[i], envir=env)
    paletteFrame <- tkframe(top)
    colorField1 <- labelRcmdr(paletteFrame, text=rgb2col(hex[1]), fg=hex[1])
    button1 <- tkbutton(paletteFrame, text=hex[1], bg = hex[1],
        fg=convert(hex[1]),
        command=function() {
            color <- pickColor(hex[1], parent=button1)
            fg <- convert(color)
            tkconfigure(button1, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField1, text=rgb2col(color), foreground=color)
            assign("hex.1", color, envir=env)
        }
    )
    colorField2 <- labelRcmdr(paletteFrame, text=rgb2col(hex[2]), fg=hex[2])
    button2 <- tkbutton(paletteFrame, text=hex[2], bg = hex[2],
        fg=convert(hex[2]),
        command=function() {
            color <- pickColor(hex[2], parent=button2)
            fg <- convert(color)
            tkconfigure(button2, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField2, text=rgb2col(color), foreground=color)
            assign("hex.2", color, envir=env)
        }
    )
    colorField3 <- labelRcmdr(paletteFrame, text=rgb2col(hex[3]), fg=hex[3])
    button3 <- tkbutton(paletteFrame, text=hex[3], bg = hex[3],
        fg=convert(hex[3]),
        command=function() {
            color <- pickColor(hex[3], parent=button3)
            fg <- convert(color)
            tkconfigure(button3, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField3, text=rgb2col(color), foreground=color)
            assign("hex.3", color, envir=env)
        }
    )
    colorField4 <- labelRcmdr(paletteFrame, text=rgb2col(hex[4]), fg=hex[4])
    button4 <- tkbutton(paletteFrame, text=hex[4], bg = hex[4],
        fg=convert(hex[4]),
        command=function() {
            color <- pickColor(hex[4], parent=button4)
            fg <- convert(color)
            tkconfigure(button4, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField4, text=rgb2col(color), foreground=color)
            assign("hex.4", color, envir=env)
        }
    )
    colorField5 <- labelRcmdr(paletteFrame, text=rgb2col(hex[5]), fg=hex[5])
    button5 <- tkbutton(paletteFrame, text=hex[5], bg = hex[5],
        fg=convert(hex[5]),
        command=function() {
            color <- pickColor(hex[5], parent=button5)
            fg <- convert(color)
            tkconfigure(button5, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField5, text=rgb2col(color), foreground=color)
            assign("hex.5", color, envir=env)
        }
    )
    colorField6 <- labelRcmdr(paletteFrame, text=rgb2col(hex[6]), fg=hex[6])
    button6 <- tkbutton(paletteFrame, text=hex[6], bg = hex[6],
        fg=convert(hex[6]),
        command=function() {
            color <- pickColor(hex[6], parent=button6)
            fg <- convert(color)
            tkconfigure(button6, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField6, text=rgb2col(color), foreground=color)
            assign("hex.6", color, envir=env)
        }
    )
    colorField7 <- labelRcmdr(paletteFrame, text=rgb2col(hex[7]), fg=hex[7])
    button7 <- tkbutton(paletteFrame, text=hex[7], bg = hex[7],
        fg=convert(hex[7]),
        command=function() {
            color <- pickColor(hex[7], parent=button7)
            fg <- convert(color)
            tkconfigure(button7, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField7, text=rgb2col(color), foreground=color)
            assign("hex.7", color, envir=env)
        }
    )
    colorField8 <- labelRcmdr(paletteFrame, text=rgb2col(hex[8]), fg=hex[8])
    button8 <- tkbutton(paletteFrame, text=hex[8], bg = hex[8],
        fg=convert(hex[8]),
        command=function() {
            color <- pickColor(hex[8], parent=button8)
            fg <- convert(color)
            tkconfigure(button8, bg=color, fg=fg, text=toupper(color))
            tkconfigure(colorField8, text=rgb2col(color), foreground=color)
            assign("hex.8", color, envir=env)
        }
    )
    onOK <- function(){
        closeDialog(top)
        palette(c(hex.1, hex.2, hex.3, hex.4, hex.5, hex.6, hex.7, hex.8))
        Message(gettextRcmdr("Color palette reset.", type="note"))
    }
    OKCancelHelp(helpSubject="palette")
    tkgrid(button1, button2, button3, button4, button5, button6, button7, button8)
    tkgrid(colorField1, colorField2, colorField3, colorField4, colorField5, 
        colorField6, colorField7, colorField8)
    tkgrid(paletteFrame)
    tkgrid(buttonsFrame, sticky="ew")
    dialogSuffix()
}

stripChart <- function () {
    defaults <- list (initial.group = NULL, initial.response = NULL, initial.plotType = "stack",
        initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"), 
        initial.main=gettextRcmdr("<auto>"), initial.tab=0)
    dialog.values <- getDialog("stripChart", defaults)
    initializeDialog(title = gettextRcmdr("Strip Chart"), use.tabs=TRUE)
    groupBox <- variableListBox(dataTab, Factors(), title = gettextRcmdr("Factors (pick zero or more)"),
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.group, "factor"))
    responseBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.response, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Duplicate Values"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
        closeDialog()
        if (0 == length(response)) {
            errorCondition(recall = stripChart, message = gettextRcmdr("No response variable selected."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        plotType <- tclvalue(plotTypeVariable)
        
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ""
        else paste(", main=\"", main, "\"", sep = "")
        putDialog ("stripChart", list (initial.group = groups, initial.response = response,
            initial.plotType = plotType, initial.xlab=tclvalue(xlabVar),
            initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar),
            initial.tab=tab))
        method <- paste(", method=\"", plotType, "\"", sep = "")
        if (length(groups) == 0) {
            xlab <- trim.blanks(tclvalue(xlabVar))
            xlab <- if (xlab == gettextRcmdr("<auto>"))
                paste(", xlab=\"", response, "\"", sep = "")
            else paste(", xlab=\"", xlab, "\"", sep = "")
            ylab <- trim.blanks(tclvalue(ylabVar))
            ylab <- if (ylab == gettextRcmdr("<auto>"))
                ""
            else paste(", ylab=\"", ylab, "\"", sep = "")
            doItAndPrint(paste("stripchart(", .activeDataSet,
                "$", response, method, xlab, ylab, main, ")", sep = ""))
        }
        else {
            groupNames <- paste(groups, collapse = "*")
            xlab <- trim.blanks(tclvalue(xlabVar))
            xlab <- if (xlab == gettextRcmdr("<auto>"))
                ""
            else paste(", xlab=\"", xlab, "\"", sep = "")
            ylab <- trim.blanks(tclvalue(ylabVar))
            ylab <- if (ylab == gettextRcmdr("<auto>"))
                paste(", ylab=\"", response, "\"", sep = "")
            else paste(", ylab=\"", ylab, "\"", sep = "")
            doItAndPrint(paste("stripchart(", response, " ~ ",
                groupNames, ", vertical=TRUE", method, xlab, ylab, main, ", data=",
                .activeDataSet, ")", sep = ""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    plotFrame <- tkframe(optFrame)
    radioButtons(window = plotFrame, name = "plotType", buttons = c("stack", "jitter"),
        labels = gettextRcmdr(c("Stack", "Jitter")), title = "",
        initialValue = dialog.values$initial.plotType)
    OKCancelHelp(helpSubject = "stripchart", reset = "stripChart", apply = "stripChart")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
    tkgrid(plotTypeFrame, sticky = "w")
    tkgrid(plotFrame, sticky = "w")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

DensityPlot <- function () {
    defaults <- list(initial.x = NULL, initial.bw = gettextRcmdr("<auto>"),
        initial.kernel="gaussian", initial.adjust=1, initial.group=NULL, initial.tab=0,
        initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"))
    dialog.values <- getDialog("DensityPlot", defaults)
    initializeDialog(title = gettextRcmdr("Nonparametric Density Estimate"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, Numeric(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    optionsFrame <- tkframe(optionsTab)
    optFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Kernel Function"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
        font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    kerFrame <- tkframe(optFrame)
    radioButtons(kerFrame, name = "kernel", buttons = c("gaussian", "epanechnikov", "biweight"),
        labels = gettextRcmdr(c("Gaussian", "Epanechnikov", "Tukey biweight")),
        initialValue = dialog.values$initial.kernel)
    bwFrame <- tkframe(optFrame)
    bwVariable <- tclVar(dialog.values$initial.bw)
    bwField <- ttkentry(bwFrame, width = "8", textvariable = bwVariable)
    adjustVariable <- tclVar(dialog.values$initial.adjust)
    adjustSlider <- tkscale(bwFrame, from = 0.1, to = 10, showvalue = TRUE,
        variable = adjustVariable, resolution = 0.1, orient = "horizontal")
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
        command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
        ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        kernel <- tclvalue(kernelVariable)
        adjust <- tclvalue(adjustVariable)
        bw <- tclvalue(bwVariable)
        if (length(x) == 0) {
            errorCondition(recall = DensityPlot, message = gettextRcmdr("You must select a variable"))
            return()
        }
        if (bw != gettextRcmdr("<auto>")){
            test.bw <- suppressWarnings(as.numeric(bw))
            if (is.na(test.bw) || test.bw <= 0){
                errorCondition(recall = DensityPlot,
                    message = gettextRcmdr("Bandwidth must be <auto> or a positive number"))
                return()
            }
        }
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            ""
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>"))
            ""
        else paste(", ylab=\"", ylab, "\"", sep = "")
        putDialog ("DensityPlot", list(initial.x = x, initial.bw = bw, initial.kernel=kernel,
            initial.adjust=adjust,
            initial.group=if (.groups == FALSE) NULL else .groups,
            initial.tab=tab, initial.xlab=tclvalue(xlabVar),
            initial.ylab=tclvalue(ylabVar)))
        if (bw == gettextRcmdr("<auto>")) bw  <- '"SJ"'
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        var <- paste(.activeDataSet, "$", x, sep = "")
        if (is.null(.groups) || .groups == FALSE) {
            command <- paste("densityPlot( ~ ", x, ", data=", .activeDataSet, ', bw=', bw,
                ", adjust=", adjust, ', kernel="', kernel, '"', xlab, ylab, ')', sep="")
            doItAndPrint(command)
        }
        else {
            command <- paste("densityPlot(", x, "~", .groups, ", data=",
                .activeDataSet, ', bw=', bw,
                ", adjust=", adjust, ', kernel="', kernel, '"', xlab, ylab, ')', sep="")
            doItAndPrint(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(DensityPlot, initialGroup=initial.group,
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups")
        else paste(gettextRcmdr("Plot by:"), initial.group),
        window=dataTab)
    OKCancelHelp(helpSubject = "densityPlot", reset = "DensityPlot", apply="DensityPlot")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(kernelFrame, stick = "w")
    tkgrid(kerFrame, sticky="w")
    tkgrid(labelRcmdr(bwFrame, text = gettextRcmdr("Bandwidth")), bwField,
        labelRcmdr(bwFrame, text = gettextRcmdr("Multiply bandwidth by")),
        adjustSlider, sticky = "swe", padx=6)
    tkgrid(bwFrame, sticky="sw")
    tkgrid(optFrame, parFrame, sticky = "nswe", padx=6, pady=6)
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}
