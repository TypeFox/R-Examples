scatterPlot.HH <- function(){
    initializeDialog(title=gettextRcmdr("Scatterplot.HH"))
    .numeric <- Numeric()
    variablesFrame <- tkframe(top)
    xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("x-variable (pick one)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("y-variable (pick one)"))
    optionsParFrame <- tkframe(top)
    checkBoxes(window=optionsParFrame, frame="optionsFrame", 
		boxes=c("identify", "jitterX", "jitterY", "logX", "logY", "boxplots", "lsLine", "smoothLine", "spread"),
        initialValues=c(0, 0, 0, 0, 0, 0, 1, 0, 0), 
		labels=gettextRcmdr(c("Identify points", "Jitter x-variable", "Jitter y-variable", "Log x-axis", "Log y-axis",
        "Marginal boxplots", "Least-squares line", "Smooth line", "Show spread")), title="Options")
    sliderValue <- tclVar("50")
    slider <- tkscale(optionsFrame, from=0, to=100, showvalue=TRUE, variable=sliderValue,
        resolution=5, orient="horizontal")
    subsetBox()
    labelsFrame <- tkframe(top)
    xlabVar <- tclVar(gettextRcmdr("<auto>"))
    ylabVar <- tclVar(gettextRcmdr("<auto>"))
    xlabFrame <- tkframe(labelsFrame)
    xlabEntry <- ttkentry(xlabFrame, width="25", textvariable=xlabVar)
    xlabScroll <- ttkscrollbar(xlabFrame, orient="horizontal",
        command=function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand=function(...) tkset(xlabScroll, ...))
    tkgrid(labelRcmdr(xlabFrame, text=gettextRcmdr("x-axis label"), fg="blue"), sticky="w")
    tkgrid(xlabEntry, sticky="w")
    tkgrid(xlabScroll, sticky="ew")
    ylabFrame <- tkframe(labelsFrame)
    ylabEntry <- ttkentry(ylabFrame, width="25", textvariable=ylabVar)
    ylabScroll <- ttkscrollbar(ylabFrame, orient="horizontal",
        command=function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand=function(...) tkset(ylabScroll, ...))
    tkgrid(labelRcmdr(ylabFrame, text=gettextRcmdr("y-axis label"), fg="blue"), sticky="w")
    tkgrid(ylabEntry, sticky="w")
    tkgrid(ylabScroll, sticky="ew")
    tkgrid(xlabFrame, labelRcmdr(labelsFrame, text="     "), ylabFrame, sticky="w")
    parFrame <- tkframe(optionsParFrame)
    pchVar <- tclVar(gettextRcmdr("16"))
    pchEntry <- ttkentry(parFrame, width=25, textvariable=pchVar)
    cexValue <- tclVar("1.3")
    cex.axisValue <- tclVar("1.3")
    cex.labValue <- tclVar("1.3")
    cexSlider <- tkscale(parFrame, from=0.5, to=2.5, showvalue=TRUE, variable=cexValue,
        resolution=0.1, orient="horizontal")
    cex.axisSlider <- tkscale(parFrame, from=0.5, to=2.5, showvalue=TRUE, variable=cex.axisValue,
        resolution=0.1, orient="horizontal")
    cex.labSlider <- tkscale(parFrame, from=0.5, to=2.5, showvalue=TRUE, variable=cex.labValue,
        resolution=0.1, orient="horizontal")
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (length(x) == 0 || length(y) == 0){
            errorCondition(recall=scatterPlot.HH, message=gettextRcmdr("You must select two variables"))
            return()
            }
        if (x == y) {
            errorCondition(recall=scatterPlot.HH, message=gettextRcmdr("x and y variables must be different"))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        jitter <- if ("1" == tclvalue(jitterXVariable) && "1" == tclvalue(jitterYVariable)) ", jitter=list(x=1, y=1)"
            else if ("1" == tclvalue(jitterXVariable)) ", jitter=list(x=1)"
            else if ("1" == tclvalue(jitterYVariable)) ", jitter=list(y=1)"
            else ""
		logstring <- ""
		if ("1" == tclvalue(logXVariable)) logstring <- paste(logstring, "x", sep="")
		if ("1" == tclvalue(logYVariable)) logstring <- paste(logstring, "y", sep="")
		log <- if(logstring != "") paste(', log="', logstring, '"', sep="") else ""
		if("1" == tclvalue(identifyVariable)){
			RcmdrTkmessageBox(title="Identify Points",
					message=paste(gettextRcmdr("Use left mouse button to identify points,\n"),
						gettextRcmdr(if (MacOSXP()) "esc key to exit." else "right button to exit."), sep=""),
					icon="info", type="ok")
			idtext <- ', id.method="identify"'
		}
        else idtext <- ""
        box <- if ("1" == tclvalue(boxplotsVariable)) "'xy'" else "FALSE"
        line <- if("1" == tclvalue(lsLineVariable)) "lm" else "FALSE"
        smooth <- as.character("1" == tclvalue(smoothLineVariable))
		spread <- as.character("1" == tclvalue(spreadVariable))
        span <- as.numeric(tclvalue(sliderValue))
        subset <- tclvalue(subsetVariable)
        subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
            else paste(", subset=", subset, sep="")
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if(xlab == gettextRcmdr("<auto>")) "" else paste(', xlab="', xlab, '"', sep="")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if(ylab == gettextRcmdr("<auto>")) "" else paste(', ylab="', ylab, '"', sep="")
        cex <- as.numeric(tclvalue(cexValue))
        cex <- if(cex == 1) "" else paste(', cex=', cex, sep="")
        cex.axis <- as.numeric(tclvalue(cex.axisValue))
        cex.axis <- if(cex.axis == 1) "" else paste(', cex.axis=', cex.axis, sep="")
        cex.lab <- as.numeric(tclvalue(cex.labValue))
        cex.lab <- if(cex.lab == 1) "" else paste(', cex.lab=', cex.lab, sep="")
        pch <- gsub(" ", ",", tclvalue(pchVar))
        if ("" == pch) {
            errorCondition(recall=scatterPlot.HH, message=gettextRcmdr("No plotting characters."))
            return()
            }
        pch <- if(trim.blanks(pch) == gettextRcmdr("<auto>")) "" else paste(", pch=c(", pch, ")", sep="")
        if (.groups == FALSE) {
            doItAndPrint(paste("scatterplot(", y, "~", x, log,
                ", reg.line=", line, ", smooth=", smooth, ", spread=", spread, idtext,
                ", boxplots=", box, ", span=", span/100, jitter, xlab, ylab,
                cex, cex.axis, cex.lab, pch,
                ", data=", .activeDataSet, subset, ")", sep=""))
            }
        else {
            doItAndPrint(paste("scatterplot(", y, "~", x," | ", .groups,
                ", reg.line=", line, ", smooth=", smooth, ", spread=", spread, idtext,
                ", boxplots=", box, ", span=", span/100, jitter, xlab, ylab,
                cex, cex.axis, cex.lab, pch,
                ", by.groups=", .linesByGroup,
                ", data=", .activeDataSet, subset, ")", sep=""))
            }
        activateMenus()
        tkfocus(CommanderWindow())
        }
    groupsBox(scatterPlot, plotLinesByGroup=TRUE)
    OKCancelHelp(helpSubject="scatterplot")
    tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Span for smooth")), slider, sticky="w")
    tkgrid(labelRcmdr(parFrame, text=gettextRcmdr("Plotting Parameters"), fg="blue"), sticky="w")
    tkgrid(labelRcmdr(parFrame, text=gettextRcmdr("Plotting characters")), pchEntry, stick="w")
    tkgrid(labelRcmdr(parFrame, text=gettextRcmdr("Point size")), cexSlider, sticky="w")
    tkgrid(labelRcmdr(parFrame, text=gettextRcmdr("Axis text size")), cex.axisSlider, sticky="w")
    tkgrid(labelRcmdr(parFrame, text=gettextRcmdr("Axis-labels text size")), cex.labSlider, sticky="w")
    tkgrid(optionsFrame, parFrame, sticky="nw")
    tkgrid(optionsParFrame, sticky="w")
    tkgrid(labelsFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(groupsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=" "))
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=8, columns=2)
    }
