

oneWayAnova2 <- function () {
	
	defaults <- list(initial.group = NULL, initial.response = NULL)
	dialog.values <- getDialog("oneWayAnova2", defaults)
	initializeDialog(title = gettextRcmdr("Significativit\xE9 N*C"))
	UpdateModelNumber()
	modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), 
					sep = ""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
	groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Groups (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.group, "factor"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	optionsFrame <- tkframe(top)
	

	
onOK <- function() {
		modelValue <- trim.blanks(tclvalue(modelName))
		if (!is.valid.name(modelValue)) {
			UpdateModelNumber(-1)
			errorCondition(recall = oneWayAnova2, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
							modelValue))
			return()
		}
		if (is.element(modelValue, listAOVModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
				UpdateModelNumber(-1)
				tkdestroy(top)
				oneWayAnova2()
				return()
			}
		}
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = oneWayAnova2, message = gettextRcmdr("You must select a groups factor."))
			return()
		}
		if (length(response) == 0) {
			errorCondition(recall = oneWayAnova2, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		command <- paste(modelValue, " <- aov(", response, " ~ ", 
				group, ", data=", .activeDataSet, ")", sep = "")
		justDoIt(command)
		logger(command)
		doItAndPrint(paste("signifNC(", modelValue, ")", sep = ""))
		
		activeModel(modelValue)

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "anova", model = TRUE, reset = "oneWayAnova2")
	tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model: ")), 
			model, sticky = "w")
	tkgrid(modelFrame, sticky = "w", columnspan = 2)
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	
	tkgrid(optionsFrame, sticky = "w", columnspan = 2)
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 4, columns = 2)
}
