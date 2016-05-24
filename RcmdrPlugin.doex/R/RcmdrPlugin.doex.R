# Some Rcmdr dialogs for the doex package

# last modified: March 28 2010 by E. Hodgess

# Note: the following function (with contributions from Richard Heiberger) 
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded

.onAttach <- function(libname, pkgname){
        if (!interactive()) return()
        Rcmdr <- options()$Rcmdr
        plugins <- Rcmdr$plugins
        if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
                Rcmdr$plugins <- c(plugins, pkgname)
                options(Rcmdr=Rcmdr)
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
        }
}



	ranblockAnova <- function(){
		initializeDialog(title=gettextRcmdr("Randomized Block Design"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick two)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			modelValue <- trim.blanks(tclvalue(modelName))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=ranblockAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					ranblockAnova()
					return()
				}
			}
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) != 2){
				errorCondition(recall=ranblockAnova, message=gettextRcmdr("You must select  two factors."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=ranblockAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste(modelValue, " <- (aov(", response, " ~ ", paste(groups, collapse="+"),
					", data=", .activeDataSet, "))", sep=""))
			doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
			doItAndPrint(paste("TukeyHSD(",modelValue,")",
sep=""))
			
					justDoIt("par(mfrow=c(2,1))")
					justDoIt("old.oma <- par(oma=c(0,5,0,0))")
					logger("old.oma <- par(oma=c(0,5,0,0))")
command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[1], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					comm1 <- paste("mtext('",groups[1],
"')",sep="")
					
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
					justDoIt(comm1)
					logger(comm1)
					logger("plot(confint(.Pairs))")
	command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[2], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					
				comm1 <- paste("mtext('",groups[2],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)
					justDoIt("par(old.oma)")
					logger("par(old.oma)")
					logger("remove(.Pairs)")
					remove(.Pairs, envir=.GlobalEnv)
					justDoIt("par(mfrow=c(1,1))")

			activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="Anova", model=TRUE)
		tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}
	



	LatinAnova <- function(){
		initializeDialog(title=gettextRcmdr("Latin Square Design"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick three)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			modelValue <- trim.blanks(tclvalue(modelName))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=LatinAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					ranblockAnova()
					return()
				}
			}
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) != 3){
				errorCondition(recall=LatinAnova, message=gettextRcmdr("You must select  three factors."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=LatinAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste(modelValue, " <- (aov(", response, " ~ ", paste(groups, collapse="+"),
					", data=", .activeDataSet, "))", sep=""))
			doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
			doItAndPrint(paste("TukeyHSD(",modelValue,")",
sep=""))
					justDoIt("par(mfrow=c(3,1))")
					justDoIt("old.oma <- par(oma=c(0,5,0,0))")
					logger("old.oma <- par(oma=c(0,5,0,0))")
command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[1], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					comm1 <- paste("mtext('",groups[1],
"')",sep="")
					
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
					justDoIt(comm1)
					logger(comm1)
					logger("plot(confint(.Pairs))")
	command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[2], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					
				comm1 <- paste("mtext('",groups[2],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)
										logger("remove(.Pairs)")
					remove(.Pairs, envir=.GlobalEnv)
					
	command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[3], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					
				comm1 <- paste("mtext('",groups[3],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)
					justDoIt("par(old.oma)")
					logger("par(old.oma)")
					logger("remove(.Pairs)")
					remove(.Pairs, envir=.GlobalEnv)
					justDoIt("par(mfrow=c(1,1))")
					justDoIt("par(mfrow=c(1,1))")


			activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="Anova", model=TRUE)
		tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}





	bartlett <- function(){
		initializeDialog(title=gettextRcmdr("Bartlett Test"))
		UpdateModelNumber()
		
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Treatment Factor (pick one)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) != 1){
				errorCondition(recall=bartlett, message=gettextRcmdr("You must select  one factor."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=ranblockAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste("bartlett.test(", response, " ~ ", groups,
					", data=", .activeDataSet, ")", sep=""))
			

	
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="bartlett.test", model=TRUE)
		
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}
	





	GreekAnova <- function(){
		initializeDialog(title=gettextRcmdr("Graeco-Latin Square Design"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick four)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			modelValue <- trim.blanks(tclvalue(modelName))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=GreekAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					ranblockAnova()
					return()
				}
			}
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) != 4){
				errorCondition(recall=GreekAnova, message=gettextRcmdr("You must select  three factors."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=GreekAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste(modelValue, " <- (aov(", response, " ~ ", paste(groups, collapse="+"),
					", data=", .activeDataSet, "))", sep=""))
			doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
			doItAndPrint(paste("TukeyHSD(",modelValue,")",
sep=""))
					justDoIt("par(mfrow=c(4,1))")
					justDoIt("old.oma <- par(oma=c(0,5,0,0))")
					logger("old.oma <- par(oma=c(0,5,0,0))")
command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[1], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					comm1 <- paste("mtext('",groups[1],
"')",sep="")
					
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
					justDoIt(comm1)
					logger(comm1)
					logger("plot(confint(.Pairs))")
	command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[2], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					
				comm1 <- paste("mtext('",groups[2],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)
										logger("remove(.Pairs)")
					remove(.Pairs, envir=.GlobalEnv)
					
	command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", groups[3], ' = "Tukey"))', sep="")
					justDoIt(command)
					logger(command)
					doItAndPrint("confint(.Pairs) # confidence intervals")
					
				comm1 <- paste("mtext('",groups[3],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)

					comm1 <- paste("mtext('",groups[4],
"')",sep="")
					justDoIt("plot(confint(.Pairs))")
					justDoIt("abline(v=0,col='red')")
										logger("plot(confint(.Pairs))")
					justDoIt(comm1)
					justDoIt("par(old.oma)")
					logger("par(old.oma)")
					logger("remove(.Pairs)")
					remove(.Pairs, envir=.GlobalEnv)
					justDoIt("par(mfrow=c(1,1))")
					justDoIt("par(mfrow=c(1,1))")


			activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="Anova", model=TRUE)
		tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}


	factorialAnova <- function(){
		initializeDialog(title=gettextRcmdr("Two Factor Factorial Design"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick two)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			modelValue <- trim.blanks(tclvalue(modelName))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=factorialAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					factorialAnova()
					return()
				}
			}
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) != 2){
				errorCondition(recall=factorialAnova, message=gettextRcmdr("You must select  two factors."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=factorialAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste(modelValue, " <- (aov(", response, " ~ ", paste(groups, collapse="*"),
					", data=", .activeDataSet, "))", sep=""))
			doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
			
			
	             doItAndPrint(paste("plot(density(residuals(",modelValue,")))",sep=""))
					
			activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="Anova", model=TRUE)
		tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}
	



	twokfactorialAnova <- function(){
		initializeDialog(title=gettextRcmdr("Two k Factor Factorial Design"))
		UpdateModelNumber()
		modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
		modelFrame <- tkframe(top)
		model <- ttkentry(modelFrame, width="20", textvariable=modelName)
		groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick at least two)"))
		responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
		onOK <- function(){
			modelValue <- trim.blanks(tclvalue(modelName))
			if (!is.valid.name(modelValue)){
				UpdateModelNumber(-1)
				errorCondition(recall=twokfactorialAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
				return()
			}
			if (is.element(modelValue, listAOVModels())) {
				if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
					UpdateModelNumber(-1)
					tkdestroy(top)
					twokfactorialAnova()
					return()
				}
			}
			groups <- getSelection(groupBox)
			response <- getSelection(responseBox)
			closeDialog()
			if (length(groups) < 2){
				errorCondition(recall=twokfactorialAnova, message=gettextRcmdr("You must select at least two factors."))
				return()
			}
			if (length(response) == 0){
				errorCondition(recall=twokfactorialAnova, message=gettextRcmdr("You must select a response variable."))
				return()
			}
			.activeDataSet <- ActiveDataSet()
			groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
			doItAndPrint(paste(modelValue, " <- (aov(", response, " ~ ", paste(groups, collapse="*"),
					", data=", .activeDataSet, "))", sep=""))
			doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
			
			
	             doItAndPrint(paste("plot(density(residuals(",modelValue,")))",sep=""))
					
			activeModel(modelValue)
			tkfocus(CommanderWindow())
		}
		OKCancelHelp(helpSubject="Anova", model=TRUE)
		tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
		tkgrid(modelFrame, sticky="w", columnspan=2)
		tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
		tkgrid(buttonsFrame, columnspan=2, sticky="w")
		dialogSuffix(rows=4, columns=2)
	}
	


