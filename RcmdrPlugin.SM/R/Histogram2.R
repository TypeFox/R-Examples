Histogram2 <- function () {
	defaults <- list(initial.x = NULL, initial.scale = "frequency", 
			initial.bins = gettextRcmdr ("<auto>")) 
	dialog.values <- getDialog("Histogram2", defaults)
	initializeDialog(title = gettextRcmdr("Histogram"))
	xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	onOK <- function() {
		x <- getSelection(xBox)
		closeDialog()
		if (length(x) == 0) {
			errorCondition(recall = Histogram2, message = gettextRcmdr("You must select a variable"))
			return()
		}
		bins <- tclvalue(binsVariable)
		opts <- options(warn = -1)
		binstext <- if (bins == gettextRcmdr("<auto>")) 
					"\"Sturges\""
				else as.numeric(bins)
		options(opts)
scale <- tclvalue(scaleVariable)

		putDialog ("Histogram2", list (initial.x = x, initial.bins = bins, initial.scale = scale))


ttt<-c("frequency","percent","density")
scale.n<-charmatch(scale,ttt)
name.scale<-c("effectif","pourcentage","densit\xE9")[scale.n]


		command <- paste("Hist(", ActiveDataSet(), "$", x, ", scale=\"", 
				scale, "\", breaks=", binstext, ", col=\"pink\",xlab=\"",x,"\",ylab=\"",name.scale,"\")", 
				sep = "")
		doItAndPrint(command)
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "Hist", reset = "Histogram2")
	radioButtons(name = "scale", buttons = c("frequency", "percent", 
					"density"), labels = gettextRcmdr(c("Frequency counts", 
							"Percentages", "Densities")), title = gettextRcmdr("Axis Scaling"), 
			initialValue = dialog.values$initial.scale)
	binsFrame <- tkframe(top)
	binsVariable <- tclVar(dialog.values$initial.bins)
	binsField <- ttkentry(binsFrame, width = "8", textvariable = binsVariable)
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(labelRcmdr(binsFrame, text = gettextRcmdr("Number of bins: ")), 
			binsField, sticky = "w")
	tkgrid(binsFrame, sticky = "w")
	tkgrid(scaleFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	tkgrid.configure(binsField, sticky = "e")
	dialogSuffix(rows = 4, columns = 1)
}


