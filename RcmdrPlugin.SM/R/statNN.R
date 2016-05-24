

statNN <- function () {
	
	defaults <- list(initial.group = NULL, initial.response = NULL)
	dialog.values <- getDialog("statNN", defaults)
	initializeDialog(title = gettextRcmdr("Statistiques N*N"))
	UpdateModelNumber()
	modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), 
					sep = ""))
	modelFrame <- tkframe(top)
	model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
	groupBox <- variableListBox(top, Numeric(), title = gettextRcmdr("x-variable (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.group, "numeric"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("y-variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	optionsFrame <- tkframe(top)
	

	
onOK <- function() {
		modelValue <- trim.blanks(tclvalue(modelName))
		if (!is.valid.name(modelValue)) {
			UpdateModelNumber(-1)
			errorCondition(recall = statNN, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
							modelValue))
			return()
		}
		if (is.element(modelValue, listAOVModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
				UpdateModelNumber(-1)
				tkdestroy(top)
				LM2()
				return()
			}
		}
		group <- getSelection(groupBox)
		response <- getSelection(responseBox)
		closeDialog()
		if (length(group) == 0) {
			errorCondition(recall = statNN, message = gettextRcmdr("You must select an explanatory variable."))
			return()
		}
		if (length(response) == 0) {
			errorCondition(recall = statNN, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		command <- paste(modelValue, " <- lm(", response, " ~ ", 
				group, ", data=", .activeDataSet, ")", sep = "")
		justDoIt(command)
		logger(command)

		doItAndPrint(paste("cor(",.activeDataSet,"$",group,",",.activeDataSet,"$",response,",use=\"complete.obs\")   # r", sep = ""))


		doItAndPrint(paste("coef(", modelValue, ")", sep = ""))
		
		activeModel(modelValue)

		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "lm", model = TRUE, reset = "statNN")
	tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model: ")), 
			model, sticky = "w")
	tkgrid(modelFrame, sticky = "w", columnspan = 2)
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	
	tkgrid(optionsFrame, sticky = "w", columnspan = 2)
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 4, columns = 2)
}
