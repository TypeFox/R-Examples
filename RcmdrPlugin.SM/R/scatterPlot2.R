
scatterPlot2 <- function () {
	require("car")
	defaults <- list(initial.x = NULL, initial.y = NULL, initial.jitterx = 0, initial.jittery = 0, 
			initial.logstringx = 0, initial.logstringy = 0, initial.log = 0, initial.box = 0, 
			initial.line = 1, initial.smooth = 1, initial.spread = 0, initial.span = 50,
			initial.subset = gettextRcmdr ("<all valid cases>"), initial.ylab = gettextRcmdr ("<auto>"), 
			initial.xlab = gettextRcmdr("<auto>"), initial.pch = gettextRcmdr("<auto>"), 
			initial.cexValue = 1, initial.cex.axisValue = 1, initial.cex.labValue = 1, initialGroup=NULL, initial.lines.by.group=1) 
	dialog.values <- getDialog("scatterPlot2", defaults)
	initial.group <- dialog.values$initial.group
	.linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	initializeDialog(title = gettextRcmdr("Scatterplot"))
	.numeric <- Numeric()
	variablesFrame <- tkframe(top)
	xBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("x-variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	yBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("y-variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.y, "numeric"))
	optionsParFrame <- tkframe(top)
	checkBoxes(window = optionsParFrame, frame = "optionsFrame", 
			boxes = c("identify", "jitterX", "jitterY", "logX", "logY", 
					"boxplots", "lsLine", "smoothLine", "spread"), initialValues = c(dialog.values$initial.log, 
					dialog.values$initial.jitterx, dialog.values$initial.jittery, 
					dialog.values$initial.logstringx, dialog.values$initial.logstringy,
					dialog.values$initial.box, dialog.values$initial.line, dialog.values$initial.smooth,
					dialog.values$initial.spread),labels = gettextRcmdr(c("Identify points", 
							"Jitter x-variable", "Jitter y-variable", "Log x-axis", 
							"Log y-axis", "Marginal boxplots", "Least-squares line", 
							"Smooth line", "Show spread")), title = gettextRcmdr("Options"))
	sliderValue <- tclVar(dialog.values$initial.span)
	slider <- tkscale(optionsFrame, from = 0, to = 100, showvalue = TRUE, 
			variable = sliderValue, resolution = 5, orient = "horizontal")
	subsetBox(subset.expression = dialog.values$initial.subset)
	labelsFrame <- tkframe(top)
	xlabVar <- tclVar(dialog.values$initial.xlab)
	ylabVar <- tclVar(dialog.values$initial.ylab)
	xlabFrame <- tkframe(labelsFrame)
	xlabEntry <- ttkentry(xlabFrame, width = "25", textvariable = xlabVar)
	xlabScroll <- ttkscrollbar(xlabFrame, orient = "horizontal", 
			command = function(...) tkxview(xlabEntry, ...))
	tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll, 
						...))
	tkgrid(labelRcmdr(xlabFrame, text = gettextRcmdr("x-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(xlabEntry, sticky = "w")
	tkgrid(xlabScroll, sticky = "ew")
	ylabFrame <- tkframe(labelsFrame)
	ylabEntry <- ttkentry(ylabFrame, width = "25", textvariable = ylabVar)
	ylabScroll <- ttkscrollbar(ylabFrame, orient = "horizontal", 
			command = function(...) tkxview(ylabEntry, ...))
	tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll, 
						...))
	tkgrid(labelRcmdr(ylabFrame, text = gettextRcmdr("y-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(ylabEntry, sticky = "w")
	tkgrid(ylabScroll, sticky = "ew")
	tkgrid(xlabFrame, labelRcmdr(labelsFrame, text = "     "), 
			ylabFrame, sticky = "w")
	parFrame <- tkframe(optionsParFrame)
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
	onOK <- function() {
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
		log <- tclvalue(identifyVariable)
		box <- tclvalue(boxplotsVariable)
		line <- tclvalue(lsLineVariable)
		smooth <-  tclvalue(smoothLineVariable)
		spread <- tclvalue(spreadVariable)
		span <- as.numeric(tclvalue(sliderValue))
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
		pch <- gsub(" ", ",", tclvalue(pchVar))
		putDialog ("scatterPlot", list (initial.x = x, initial.y = y, initial.jitterx = tclvalue(jitterXVariable),
						initial.jittery = tclvalue(jitterYVariable), initial.logstringx = tclvalue(logXVariable),
						initial.logstringy = tclvalue(logYVariable), initial.log = log, initial.box = box, 
						initial.line = line, initial.smooth = smooth, initial.spread = spread,
						initial.span = span, initial.subset = initial.subset, initial.xlab = tclvalue(xlabVar),
						initial.ylab = tclvalue(ylabVar), initial.cexValue = tclvalue(cexValue), 
						initial.cex.axisValue = tclvalue(cex.axisValue), initial.cex.labValue = tclvalue(cex.labValue), 
						initial.pch = pch, initial.group=if (.groups == FALSE) NULL else .groups,
						initial.lines.by.group=if (.linesByGroup) 1 else 0))
		closeDialog()
		if ("" == pch) {
			errorCondition(recall = scatterPlot2, message = gettextRcmdr("No plotting characters."))
			return()
		}
		pch <- if (trim.blanks(pch) == gettextRcmdr("<auto>")) 
					""
				else paste(", pch=c(", pch, ")", sep = "")
		if (length(x) == 0 || length(y) == 0) {
			errorCondition(recall = scatterPlot2, message = gettextRcmdr("You must select two variables"))
			return()
		}
		if (x == y) {
			errorCondition(recall = scatterPlot2, message = gettextRcmdr("x and y variables must be different"))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		log <- if (logstring != "") 
					paste(", log=\"", logstring, "\"", sep = "")
				else ""
		if ("1" == tclvalue(identifyVariable)) {
			RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
							gettextRcmdr(if (MacOSXP()) 
												"esc key to exit."
											else "right button to exit."), sep = ""), icon = "info", 
					type = "ok")
			idtext <- ", id.method=\"identify\""
		}
		else idtext <- ""
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
			doItAndPrint(paste("scatterplot(", y, "~", x, log, 
							", reg.line=", line, ", smooth=", smooth, ", spread=", 
							spread, idtext, ", boxplots=", box, ", span=", 
							span/100, jitter, xlab, ylab, cex, cex.axis, 
							cex.lab, pch, ", data=", .activeDataSet, subset, 
							")", sep = ""))
		}
		else {
			doItAndPrint(paste("scatterplot(", y, "~", x, " | ", 
							.groups, log, ", reg.line=", line, ", smooth=", smooth, 
							", spread=", spread, idtext, ", boxplots=", box, 
							", span=", span/100, jitter, xlab, ylab, cex, 
							cex.axis, cex.lab, pch, ", by.groups=", .linesByGroup, 
							", data=", .activeDataSet, subset, ")", sep = ""))
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
			initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "scatterplot", reset = "scatterPlot2")
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	tkgrid(variablesFrame, sticky = "w")
	tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Span for smooth")), 
			slider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Plotting Parameters"), 
					fg = "blue"), sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Plotting characters")), 
			pchEntry, stick = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Point size")), 
			cexSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis text size")), 
			cex.axisSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis-labels text size")), 
			cex.labSlider, sticky = "w")
	tkgrid(optionsFrame, parFrame, sticky = "nw")
	tkgrid(optionsParFrame, sticky = "w")
	tkgrid(labelsFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(labelRcmdr(top, text = " "))
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 8, columns = 2)
}

