# Model menu dialogs

# last modified 2015-09-04 by J. Fox

selectActiveModel <- function(){
	models <- listAllModels()
	.activeModel <- ActiveModel()
	if ((length(models) == 1) && !is.null(.activeModel)) {
		Message(message=gettextRcmdr("There is only one model in memory."),
				type="warning")
		tkfocus(CommanderWindow())
		return()
	}
	if (length(models) == 0){
		Message(message=gettextRcmdr("There are no models from which to choose."),
				type="error")
		tkfocus(CommanderWindow())
		return()
	}
	initializeDialog(title=gettextRcmdr("Select Model"))
	.activeDataSet <- ActiveDataSet()
	initial <- if (is.null(.activeModel)) NULL else which(.activeModel == models) - 1
	modelsBox <- variableListBox(top, models, title=gettextRcmdr("Models (pick one)"), 
			initialSelection=initial)
	onOK <- function(){
		model <- getSelection(modelsBox)
		closeDialog()
		if (length(model) == 0) {
			tkfocus(CommanderWindow())
			return()
		}
		dataSet <- as.character(get(model)$call$data)
		if (length(dataSet) == 0){
			errorCondition(message=gettextRcmdr("There is no dataset associated with this model."))
			return()
		}
		dataSets <- listDataSets()
		if (!is.element(dataSet, dataSets)){
			errorCondition(message=sprintf(gettextRcmdr("The dataset associated with this model, %s, is not in memory."), dataSet))
			return()
		}
		if (is.null(.activeDataSet) || (dataSet != .activeDataSet)) activeDataSet(dataSet)
		putRcmdr("modelWithSubset", "subset" %in% names(get(model)$call))
		activeModel(model)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	nameFrame <- tkframe(top)
	tkgrid(labelRcmdr(nameFrame, fg=getRcmdr("title.color"), font="RcmdrTitleFont", text=gettextRcmdr("Current Model: ")), 
			labelRcmdr(nameFrame, text=tclvalue(getRcmdr("modelName"))), sticky="w")
	tkgrid(nameFrame, sticky="w", columnspan="2")
	tkgrid(getFrame(modelsBox), columnspan="2", sticky="w")
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	dialogSuffix()
}

# summarizeModel <- function(){
# 	.activeModel <- ActiveModel()
# 	if (is.null(.activeModel) || !checkMethod("summary", .activeModel)) return()
# 	doItAndPrint(paste("summary(", .activeModel, ", cor=FALSE)", sep=""))
# }

summarizeModel <- function(){
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("summary", .activeModel)) return()
    if (!lmP()) doItAndPrint(paste("summary(", .activeModel, ", cor=FALSE)", sep=""))
    else if (!packageAvailable("car") || !packageAvailable("sandwich")){
        doItAndPrint(paste("summary(", .activeModel, ")", sep=""))
    }
    else{
        Library("sandwich")
        defaults <- list(initial.sandwich="0", initial.type="HC3")
        dialog.values <- getDialog("summarizeModel", defaults)
        initializeDialog(title = gettextRcmdr("Linear Model Summary"))
        sandwichVar <- tclVar(dialog.values$initial.sandwich)
        sandwichFrame <- tkframe(top)
        sandwichCheckFrame <- tkframe(sandwichFrame)
        sandwichCheckBox <- ttkcheckbutton(sandwichCheckFrame, variable = sandwichVar)
        radioButtons(sandwichFrame, name = "type", buttons = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
            labels = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
            title = gettextRcmdr("Sandwich estimator"), initialValue = dialog.values$initial.type)
        onOK <- function() {
            sandwich <- tclvalue(sandwichVar)
            type <- tclvalue(typeVariable)
            putDialog ("summarizeModel", list(initial.sandwich=sandwich, initial.type = type))
            closeDialog()
            command <- if (sandwich == "1") 
                paste("summarySandwich(", .activeModel, ', type="', tolower(type), '")', sep="")
            else paste("summary(", .activeModel, ")", sep="")
            doItAndPrint(command)
        }
        OKCancelHelp(helpSubject = "summarySandwich", reset = "summarizeModel", apply="summarizeModel")
        tkgrid(sandwichCheckBox, labelRcmdr(sandwichFrame, 
            text=gettextRcmdr("Use sandwich estimator of    \ncoefficient standard errors")), 
            sticky="nw")
        tkgrid(sandwichCheckFrame, typeFrame, sticky="nw")
        tkgrid(sandwichFrame, sticky = "w")
        tkgrid(buttonsFrame, sticky = "w")
        dialogSuffix()
    }
}

plotModel <- function(){
	.activeModel <- ActiveModel()
	if (is.null(.activeModel) || !checkMethod("plot", .activeModel)) return()
	command <- "oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))"
	justDoIt(command)
	logger(command)
	doItAndPrint(paste("plot(", .activeModel, ")", sep=""))
	command <- "par(oldpar)"
	justDoIt(command)
	logger(command)
}

CRPlots <- function(){
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("crPlot", .activeModel)) return()
    
    defaults <- list(initial.span=50)
    dialog.values <- getDialog("CRPlots", defaults)
    initializeDialog(title = gettextRcmdr("Component+Residual Plots"))
    sliderValue <- tclVar(dialog.values$initial.span)
    sliderFrame <- tkframe(top)
    slider <- tkscale(sliderFrame, from = 5, to = 100, showvalue = TRUE, 
                      variable = sliderValue, resolution = 5, orient = "horizontal")
    onOK <- function(){
        span <- as.numeric(tclvalue(sliderValue))
        closeDialog()
        putDialog ("CRPlots", list(initial.span=span))
        doItAndPrint(paste("crPlots(", .activeModel, ", span=", span/100, ")", sep=""))
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "crPlots", reset = "CRPlots", apply = "CRPlots")
    tkgrid(labelRcmdr(sliderFrame, text=gettextRcmdr("Span for smooth")), slider, sticky="sw")
    tkgrid(sliderFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}

AVPlots <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("avPlot", .activeModel)) 
        return()
    defaults <- list (initial.identify = "auto", initial.id.n="2")
    dialog.values <- getDialog ("AVPlots", defaults)
    initializeDialog(title = gettextRcmdr("Added-Variable Plots"))
    identifyPointsFrame <- tkframe(top)
    radioButtons(identifyPointsFrame, name = "identify", buttons = c("auto", "mouse", 
                                                                     "not"), labels = gettextRcmdr(c("Automatically", 
                                                                                                     "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"), 
                 initialValue = dialog.values$initial.identify)    
    id.n.Var <- tclVar(dialog.values$initial.id.n) 
    npointsSpinner <- tkspinbox(identifyPointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)      
    onOK <- function() {
        id.n <- tclvalue(id.n.Var)
        identify <- tclvalue(identifyVariable)
        method <- if (identify == "mouse") "identify" else "mahal"
        id.n.use <- if (identify == "not") 0 else id.n   
        closeDialog()
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = AVPlots, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        putDialog ("AVPlots", list (initial.identify = identify, initial.id.n=id.n))
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
                                                                         gettextRcmdr(if (MacOSXP()) 
                                                                             "esc key to exit."
                                                                                      else "right button to exit."), sep = ""), icon = "info", 
                              type = "ok")
        }
        command <- paste("avPlots(", .activeModel, ', id.method="', method, '", id.n=', id.n.use,  ")", sep = "")
        if (identify == "mouse") command <- suppressMarkdown(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "avPlots", reset = "AVPlots", apply = "AVPlots")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(identifyPointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(identifyPointsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

InfluencePlot <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("influencePlot", .activeModel)) 
        return()
    defaults <- list (initial.identify = "auto", initial.id.n="2")
    dialog.values <- getDialog ("InfluencePlot", defaults)
    initializeDialog(title = gettextRcmdr("Influence Plot"))
    identifyPointsFrame <- tkframe(top)
    radioButtons(identifyPointsFrame, name = "identify", buttons = c("auto", "mouse"), labels = gettextRcmdr(c("Automatically", 
                                                                                                               "Interactively with mouse")), title = gettextRcmdr("Identify Points"), 
                 initialValue = dialog.values$initial.identify)    
    id.n.Var <- tclVar(dialog.values$initial.id.n) 
    npointsSpinner <- tkspinbox(identifyPointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)      
    onOK <- function() {
        id.n <- tclvalue(id.n.Var)
        identify <- tclvalue(identifyVariable)
        method <- if (identify == "mouse") "identify" else "noteworthy"
        closeDialog()
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = InfluencePlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        putDialog ("InfluencePlot", list (initial.identify = identify, initial.id.n=id.n))
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
                                                                         gettextRcmdr(if (MacOSXP()) 
                                                                             "esc key to exit."
                                                                                      else "right button to exit."), sep = ""), icon = "info", 
                              type = "ok")
        }
        command <- paste("influencePlot(", .activeModel, ', id.method="', method, '", id.n=', id.n,  ")", sep = "")
        if (identify == "mouse") command <- suppressMarkdown(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "influencePlot", reset = "InfluencePlot", apply = "InfluencePlot")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(identifyPointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(identifyPointsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

anovaTable <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) 
        return()
    defaults <- list (initial.type = "II", initial.sandwich="0", initial.sandwich.type="HC3")
    dialog.values <- getDialog ("anovaTable", defaults)
    initializeDialog(title = gettextRcmdr("ANOVA Table"))
    radioButtons(name = "type", buttons = c("I", "II", "III"), 
        values = c("I", "II", "III"), labels = gettextRcmdr(c("Sequential (\"Type I\")", 
            "Partial obeying marginality (\"Type II\")", "Partial ignoring marginality (\"Type III\")")), 
        title = gettextRcmdr("Type of Tests"), initialValue = dialog.values$initial.type)
    sandwichVar <- tclVar(dialog.values$initial.sandwich)
    sandwichFrame <- tkframe(top)
    sandwichCheckFrame <- tkframe(sandwichFrame)
    sandwichCheckBox <- ttkcheckbutton(sandwichCheckFrame, variable = sandwichVar)
    radioButtons(sandwichFrame, name = "sandwichType", buttons = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
        labels = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
        title = gettextRcmdr("Sandwich estimator"), initialValue = dialog.values$initial.sandwich.type)
    onOK <- function() {
        type <- as.character(tclvalue(typeVariable))
        sandwich <- tclvalue(sandwichVar)
        sandwich.type <- tclvalue(sandwichTypeVariable)
        putDialog ("anovaTable", list(initial.type = type, initial.sandwich=sandwich, 
            initial.sandwich.type=sandwich.type))
        closeDialog()
        if (is.glm <- glmP()) {
            family <- eval(parse(text = paste(.activeModel, "$family$family", 
                sep = "")))
        }
        if (type == "I") {
            if (!checkMethod("anova", .activeModel)) {
                errorCondition(message = gettextRcmdr("There is no appropriate anova method for a model of this class."))
                return()
            }
            if (sandwich == "1" && lmP()) {
                errorCondition(recall = anovaTable, 
                    message = gettextRcmdr("sandwich covariance matrix unavailable with type I tests"))
                return()
            }
            if (is.glm) {
                test <- if (family %in% c("binomial", "poisson")) 
                    "Chisq"
                else "F"
                doItAndPrint(paste("anova(", .activeModel, ", test=\"", 
                    test, "\")", sep = ""))
            }
            else doItAndPrint(paste("anova(", .activeModel, ")", 
                sep = ""))
        }
        else {
            if (!checkMethod("Anova", .activeModel)) {
                errorCondition(message = gettextRcmdr("There is no appropriate Anova method for a model of this class."))
                return()
            }
            if (is.glm) {
                test <- if (family %in% c("binomial", "poisson")) 
                    "LR"
                else "F"
                doItAndPrint(paste("Anova(", .activeModel, ", type=\"", 
                    type, "\", test=\"", test, "\")", sep = ""))
            }
            else if (lmP()){
                vcov <- if (sandwich == "1"){
                    if (sandwich.type == "HAC") paste(", vcov=vcovHAC(", .activeModel, ")", sep="")
                    else paste(", vcov=hccm(", .activeModel, ', type="', tolower(sandwich.type), '")', sep="")
                }
                else ""
                doItAndPrint(paste("Anova(", .activeModel, ", type=\"", 
                    type, "\"", vcov, ")", sep = ""))
            }
            else doItAndPrint(paste("Anova(", .activeModel, ", type=\"", 
                type, "\")", sep = ""))
            if (type == "III") 
                Message(message = gettextRcmdr("Type III tests require careful attention to contrast coding."), 
                    type = "warning")
        }
    }
    OKCancelHelp(helpSubject = "Anova", reset = "anovaTable")
    tkgrid(typeFrame, sticky = "w")
    if (lmP()){
        tkgrid(labelRcmdr(sandwichFrame, text=""))
        tkgrid(sandwichCheckBox, labelRcmdr(sandwichFrame, 
            text=gettextRcmdr("Use sandwich estimator of\ncoefficient covariance matrix   ")), 
            sticky="nw")
        tkgrid(sandwichCheckFrame, sandwichTypeFrame, sticky="nw")
        tkgrid(sandwichFrame, sticky = "w")
    }
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

VIF <- function(){
	.activeModel <- ActiveModel()
	if (is.null(.activeModel) || !checkMethod("vif", .activeModel, default=TRUE)) return()
	doItAndPrint(paste("vif(", .activeModel, ")", sep=""))
}

addObservationStatistics <- function () {
  .activeDataSet <- ActiveDataSet()
  .activeModel <- ActiveModel()
  if (is.null(.activeModel)) 
    return()
  addVariable <- function(name) {
    variable <- paste(name, ".", .activeModel, sep = "")
    if (is.element(variable, .variables)) {
      ans <- checkReplace(variable)
      if (tclvalue(ans) == "no") 
        return()
    }
    paste(variable, " <- ", name, "(", .activeModel, ")", sep = "")
  }
  if (getRcmdr("modelWithSubset")) {
    Message(message = gettextRcmdr("Observation statistics not available\nfor a model fit to a subset of the data."), 
            type = "error")
    tkfocus(CommanderWindow())
    return()
  }
  defaults <- list (initial.fitted = 1, initial.residuals = 1, initial.rstudent = 1, 
                    initial.hatvalues = 1, initial.cookd = 1, initial.obsNumbers = 1)
  dialog.values <- getDialog ("addObservationStatistics", defaults)
  initializeDialog(title = gettextRcmdr("Add Observation Statistics to Data"))
  .variables <- Variables()
  obsNumberExists <- is.element("obsNumber", .variables)
  activate <- c(checkMethod("fitted", .activeModel, default = TRUE, 
                            reportError = FALSE), checkMethod("residuals", .activeModel, 
                                                              default = TRUE, reportError = FALSE), checkMethod("rstudent", 
                                                                                                                .activeModel, reportError = FALSE), checkMethod("hatvalues", 
                                                                                                                                                                .activeModel, reportError = FALSE), checkMethod("cooks.distance", 
                                                                                                                                                                                                                .activeModel, reportError = FALSE))
  checkBoxes(frame = "selectFrame", boxes = c(c("fitted", "residuals", 
                                                "rstudent", "hatvalues", "cookd")[activate], "obsNumbers"), 
             labels = c(gettextRcmdr(c("Fitted values", "Residuals", 
                                       "Studentized residuals", "Hat-values", "Cook's distances"))[activate], 
                        gettextRcmdr("Observation indices")), initialValues = c(dialog.values$initial.fitted, 
                                                                                dialog.values$initial.residuals, dialog.values$initial.rstudent, 
                                                                                dialog.values$initial.hatvalues, dialog.values$initial.cookd, dialog.values$initial.obsNumbers))
  command <- paste(.activeDataSet, "<- within(", .activeDataSet, ", {", sep="")
  onOK <- function() {
    closeDialog()
    if (activate[1] && tclvalue(fittedVariable) == 1) 
      command <- paste(command, "\n  ", addVariable("fitted"), sep="")
    if (activate[2] && tclvalue(residualsVariable) == 1) 
      command <- paste(command, "\n  ", addVariable("residuals"), sep="")
    if (activate[3] && tclvalue(rstudentVariable) == 1) 
      command <- paste(command, "\n  ", addVariable("rstudent"), sep="")
    if (activate[4] && tclvalue(hatvaluesVariable) == 1) 
      command <- paste(command, "\n  ", addVariable("hatvalues"), sep="")
    if (activate[5] && tclvalue(cookdVariable) == 1) 
      command <- paste(command, "\n  ", addVariable("cooks.distance"), sep="")
    obsNumbers <- tclvalue(obsNumbersVariable)
    putDialog ("addObservationStatistics", list (initial.fitted = tclvalue (fittedVariable),
                                                 initial.residuals = tclvalue (residualsVariable), initial.rstudent = tclvalue(rstudentVariable), 
                                                 initial.hatvalues = tclvalue (hatvaluesVariable), initial.cookd = tclvalue (cookdVariable), 
                                                 initial.obsNumbers = obsNumbers))
    if (tclvalue(obsNumbersVariable) == 1) {
      proceed <- if (obsNumberExists) 
        tclvalue(checkReplace("obsNumber"))
      else "yes"  
      if (proceed == "yes") {
        command <- paste(command, "\n  obsNumber <- 1:nrow(", .activeDataSet, ")", sep = "")
      }
    }
    command <- paste(command, "\n})")
    result <- doItAndPrint(command)
    if (class(result) != "try-error")activeDataSet(.activeDataSet, flushModel = FALSE, flushDialogMemory = FALSE)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "influence.measures", reset = "addObservationStatistics")
  tkgrid(selectFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

residualQQPlot <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("qqPlot", .activeModel)) 
        return()
    defaults <- list (initial.simulate = 1, initial.identify = "auto", initial.id.n="2")
    dialog.values <- getDialog ("residualQQPlot", defaults)
    initializeDialog(title = gettextRcmdr("Residual Quantile-Comparison Plot"))
    selectFrame <- tkframe(top)
    simulateVar <- tclVar(dialog.values$initial.simulate)
    simulateCheckBox <- ttkcheckbutton(selectFrame, variable = simulateVar)
    identifyPointsFrame <- tkframe(top)
    radioButtons(identifyPointsFrame, name = "identify", buttons = c("auto", "mouse", 
                                                                     "not"), labels = gettextRcmdr(c("Automatically", 
                                                                                                     "Interactively with mouse", "Do not identify")), title = gettextRcmdr("Identify Points"), 
                 initialValue = dialog.values$initial.identify)    
    id.n.Var <- tclVar(dialog.values$initial.id.n) 
    npointsSpinner <- tkspinbox(identifyPointsFrame, from=1, to=10, width=2, textvariable=id.n.Var)      
    onOK <- function() {
        simulate <- tclvalue (simulateVar)  
        id.n <- tclvalue(id.n.Var)
        identify <- tclvalue(identifyVariable)
        method <- if (identify == "mouse") "identify" else "y"
        id.n.use <- if (identify == "not") 0 else id.n   
        closeDialog()
        if (is.na(suppressWarnings(as.numeric(id.n))) || round(as.numeric(id.n)) != as.numeric(id.n)){
            errorCondition(recall = residualQQPlot, message = gettextRcmdr("number of points to identify must be an integer"))
            return()
        }
        putDialog ("residualQQPlot", list (initial.simulate = simulate, initial.identify = identify, initial.id.n=id.n))
        simulate <- tclvalue(simulateVar) == 1
        if (identify == "mouse") {
            RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
                                                                         gettextRcmdr(if (MacOSXP()) 
                                                                             "esc key to exit."
                                                                                      else "right button to exit."), sep = ""), icon = "info", 
                              type = "ok")
        }
        command <- paste("qqPlot(", .activeModel, ", simulate=", 
                         simulate, ', id.method="', method, '", id.n=', id.n.use,  ")", sep = "")
        if (identify == "mouse") command <- suppressMarkdown(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qqPlot.lm", reset = "residualQQPlot", apply = "residualQQPlot")
    tkgrid(labelRcmdr(selectFrame, text = gettextRcmdr("Simulated confidence envelope")), 
           simulateCheckBox, sticky = "w")
    tkgrid(selectFrame, sticky = "w")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(labelRcmdr(identifyPointsFrame, text=gettextRcmdr("Number of points to identify  ")), npointsSpinner, sticky="w")
    tkgrid(identifyPointsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

testLinearHypothesis <- function(){
    defaults <- list(previous.model=NULL, nrows=1, table.values=0, rhs.values=0, initial.sandwich="0", initial.sandwich.type="HC3")
    dialog.values <- getDialog("testLinearHypothesis", defaults=defaults)
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("linearHypothesis", .activeModel, default=TRUE)) return()
    if (!is.null(dialog.values$previous.model)){
        if (dialog.values$previous.model != .activeModel){
            dialog.values <- defaults
        }
    }
    table.values <- dialog.values$table.values
    rhs.values <- dialog.values$rhs.values
    env <- environment()
    initializeDialog(title=gettextRcmdr("Test Linear Hypothesis"))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    setUpTable <- function(...){
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- as.numeric(tclvalue(rowsValue))
        if (length(table.values) == 1 && table.values == 0) {
            table.values <- matrix(0, nrows, ncols)
            rhs.values <- rep(0, nrows)
        }
        if (nrow(table.values) < nrows){
            add.rows <- nrows - nrow(table.values)
            table.values <- rbind(table.values, matrix(0, add.rows, ncols))
            rhs.values <- c(rhs.values, rep(0, add.rows))
        }
        col.names <- names(Coef(get(.activeModel)))
        col.names <- substring(paste(abbreviate(col.names, 12), "            "), 1, 12)
        make.col.names <- "labelRcmdr(.tableFrame, text='')"
        for (j in 1:ncols) {
            make.col.names <- paste(make.col.names, ", ", 
                                    "labelRcmdr(.tableFrame, text='", col.names[j], "')", sep="")
        }
        rhsText <- gettextRcmdr("Right-hand side")
        make.col.names <- paste(make.col.names, ", labelRcmdr(.tableFrame, text='          ')",
                                ", labelRcmdr(.tableFrame, text='", rhsText, "')", sep="")
        eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows){   
            varname <- paste(".tab.", i, ".1", sep="") 
            rhs.name <- paste(".rhs.", i, sep="")
            assign(varname, tclVar(table.values[i, 1]) , envir=env)
            assign(rhs.name, tclVar(rhs.values[i]), envir=env)
            make.row <- paste("labelRcmdr(.tableFrame, text=", i, ")")
            make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=", 
                              varname, ")", sep="")
            for (j in 2:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(table.values[i, j]), envir=env)
                make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=", 
                                  varname, ")", sep="")
            }
            make.row <- paste(make.row, ", labelRcmdr(.tableFrame, text='     '),",
                              "ttkentry(.tableFrame, width='5', textvariable=", rhs.name, ")", sep="")
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
        }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
    }
    ncols <- length(Coef(get(.activeModel)))
    rowsFrame <- tkframe(top)
    rowsValue <- tclVar(dialog.values$nrows)
    rowsSlider <- tkscale(rowsFrame, from=1, to=ncols, showvalue=FALSE, variable=rowsValue,
                          resolution=1, orient="horizontal", command=setUpTable)
    rowsShow <- labelRcmdr(rowsFrame, textvariable=rowsValue, width=2, justify="right")
    sandwichVar <- tclVar(dialog.values$initial.sandwich)
    sandwichFrame <- tkframe(top)
    sandwichCheckFrame <- tkframe(sandwichFrame)
    sandwichCheckBox <- ttkcheckbutton(sandwichCheckFrame, variable = sandwichVar)
    radioButtons(sandwichFrame, name = "sandwichType", buttons = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
                 labels = c("HC0", "HC1", "HC2", "HC3", "HC4", "HAC"), 
                 title = gettextRcmdr("Sandwich estimator"), initialValue = dialog.values$initial.sandwich.type)
    onOK <- function(){
        nrows <- as.numeric(tclvalue(rowsValue))
        cell <- 0
        values <- rep(NA, nrows*ncols)
        rhs <- rep(NA, nrows)
        for (i in 1:nrows){
            rhs.name <- paste(".rhs.", i, sep="")
            rhs[i] <- as.numeric(eval(parse(text=paste("tclvalue(", rhs.name,")", sep=""))))
            for (j in 1:ncols){
                cell <- cell+1
                varname <- paste(".tab.", i, ".", j, sep="")
                values[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
            }
        }
        values <- na.omit(values)
        sandwich <- tclvalue(sandwichVar)
        sandwich.type <- tclvalue(sandwichTypeVariable)
        closeDialog()
        if (length(values) != nrows*ncols){
            Message(message=
                        sprintf(gettextRcmdr("Number of valid entries in hypothesis matrix(%d)\nnot equal to number of rows (%d) * number of columns (%d)."), 
                                length(values), nrows, ncols), type="error")
            testLinearHypothesis()
            return()
        }
        if (qr(matrix(values, nrows, ncols, byrow=TRUE))$rank < nrows) {
            Message(message=gettextRcmdr("Hypothesis matrix is not of full row rank."),
                    type="error")
            testLinearHypothesis()
            return()
        }            
        rhs <- na.omit(rhs)
        if (length(rhs) != nrows){
            errorCondition(recall=testLinearHypothesis, 
                           message=sprintf(gettextRcmdr("Number of valid entries in rhs vector (%d)\nis not equal to number of rows (%d)."), 
                                           length(rhs), nrows))
            return()
        }
        vcov <- if (lmP() && sandwich == "1"){
            if (sandwich.type == "HAC") paste(", vcov=vcovHAC(", .activeModel, ")", sep="")
            else paste(", vcov=hccm(", .activeModel, ', type="', tolower(sandwich.type), '")', sep="")
        }
        else ""
        test <- if (glmP()) {
            family <- eval(parse(text = paste(.activeModel, "$family$family", 
                                              sep = "")))
            if (family %in% c("binomial", "poisson")) ', test="Chisq"' else ', test="F"'
        }
        else ""
        command.1 <- paste(".Hypothesis <- matrix(c(", paste(values, collapse=","), "), ", nrows, ", ", ncols,
                           ", byrow=TRUE)", sep="")
        command.2 <- paste(".RHS <- c(", paste(rhs, collapse=","), ")", sep="")
        justDoIt(paste("putRcmdr('.RHS', c(", paste(rhs, collapse=","), "))", sep=""))
        command.3 <- paste("linearHypothesis(", .activeModel, ", .Hypothesis, rhs=.RHS", vcov, test, ")", sep="")
        doItAndPrint(paste("local({\n", "  ", command.1, "\n",
                           "  ", command.2, "\n",
                           "  ", command.3, "\n",
                           "})", sep=""))                   
        tkfocus(CommanderWindow())
        contrast.table <- matrix(values, nrows, ncols, byrow=TRUE)
        putDialog("testLinearHypothesis", list(previous.model=.activeModel, nrows=nrows, table.values=contrast.table,
                                               rhs.values=getRcmdr(".RHS"), initial.sandwich=sandwich, initial.sandwich.type=sandwich.type))
    }
    OKCancelHelp(helpSubject="linearHypothesis", reset="testLinearHypothesis", apply="testLinearHypothesis")
    tkgrid(labelRcmdr(rowsFrame, text=gettextRcmdr("Number of Rows:")), rowsSlider, rowsShow, sticky="w")
    tkgrid(rowsFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter hypothesis matrix and right-hand side vector:"), 
                      fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(outerTableFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=""))
    if (lmP()){
        tkgrid(labelRcmdr(sandwichFrame, text=""))
        tkgrid(sandwichCheckBox, labelRcmdr(sandwichFrame, 
                                            text=gettextRcmdr("Use sandwich estimator of\ncoefficient covariance matrix   ")), 
               sticky="nw")
        tkgrid(sandwichCheckFrame, sandwichTypeFrame, sticky="nw")
        tkgrid(sandwichFrame, sticky = "w")
    }
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()       
} 

compareModels <- function () {
    modelPosn <- function(model){
        if (is.null(model)) return(NULL)
        if (!(model %in% models)) NULL
        else which(model == models) - 1
    }
    defaults <- list (initial.model1 = NULL, initial.model2 = NULL)
    dialog.values <- getDialog ("compareModels", defaults)  
    models <- listAllModels()
    if (length(models) < 2) {
        Message(message = gettextRcmdr("There are fewer than two models."), 
                type = "error")
        tkfocus(CommanderWindow())
        return()
    }
    initializeDialog(title = gettextRcmdr("Compare Models"))
    modelsBox1 <- variableListBox(top, models, title = gettextRcmdr("First model (pick one)"),
                                  initialSelection = modelPosn(dialog.values$initial.model1))
    modelsBox2 <- variableListBox(top, models, title = gettextRcmdr("Second model (pick one)"),
                                  initialSelection = modelPosn(dialog.values$initial.model2))
    onOK <- function() {
        model1 <- getSelection(modelsBox1)
        model2 <- getSelection(modelsBox2)
        closeDialog()
        putDialog ("compareModels", list (initial.model1 = model1, initial.model2 = model2))
        if (length(model1) == 0 || length(model2) == 0) {
            errorCondition(recall = compareModels, message = gettextRcmdr("You must select two models."))
            return()
        }
        if (!checkMethod("anova", model1)) {
            return()
        }
        if (!class(get(model1, envir = .GlobalEnv))[1] == class(get(model2, 
                                                                    envir = .GlobalEnv))[1]) {
            Message(message = gettextRcmdr("Models are not of the same class."), 
                    type = "error")
            compareModels()
            return()
        }
        if (glmP()) {
            family1 <- eval(parse(text = paste(model1, "$family$family", 
                                               sep = "")))
            family2 <- eval(parse(text = paste(model2, "$family$family", 
                                               sep = "")))
            if (family1 != family2) {
                Message(message = gettextRcmdr("Models do not have the same family."), 
                        type = "error")
                compareModels()
                return()
            }
            test <- if (family1 %in% c("binomial", "poisson")) 
                "Chisq"
            else "F"
            doItAndPrint(paste("anova(", model1, ", ", model2, 
                               ", test=\"", test, "\")", sep = ""))
        }
        else doItAndPrint(paste("anova(", model1, ", ", model2, 
                                ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "anova", reset = "compareModels", apply = "compareModels")
    tkgrid(getFrame(modelsBox1), getFrame(modelsBox2), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix()
}

BreuschPaganTest <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) 
        return()
    Library("lmtest")
    currentModel <- FALSE
    defaults <- list (initial.var = "fitted", initial.student = 0)
    dialog.values <- getDialog ("BreuschPaganTest", defaults)
    initializeDialog(title = gettextRcmdr("Breusch-Pagan Test"))
    tkgrid(labelRcmdr(top, text = gettextRcmdr("Score Test for Nonconstant Error Variance"), 
                      fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    optionsFrame <- tkframe(top)
    onOK <- function() {
        student <- tclvalue(studentVariable)
        var <- tclvalue(varVariable)
        putDialog ("BreuschPaganTest", list (initial.var = var, initial.student = student))
        type <- if (var == "fitted") 
            paste(", varformula = ~ fitted.values(", .activeModel, 
                  ")", sep = "")
        else if (var == "predictors") 
            ""
        else paste(", varformula = ~", tclvalue(rhsVariable), 
                   sep = "")
        model.formula <- as.character(formula(get(.activeModel)))
        model.formula <- paste(model.formula[2], "~", model.formula[3])
        closeDialog()
        student <- if (tclvalue(studentVariable) == 1) 
            "TRUE"
        else "FALSE"
        command <- paste("bptest(", model.formula, type, ", studentize=", 
                         student, ", data=", ActiveDataSet(), ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "bptest", reset = "BreuschPaganTest", apply = "BreuschPaganTest")
    studentVariable <- tclVar(dialog.values$initial.student)
    studentFrame <- tkframe(optionsFrame)
    studentCheckBox <- ttkcheckbutton(studentFrame, variable = studentVariable)
    tkgrid(labelRcmdr(studentFrame, text = gettextRcmdr("Studentized test statistic"), 
                      justify = "left"), studentCheckBox, sticky = "w")
    tkgrid(studentFrame, sticky = "w")
    radioButtons(optionsFrame, name = "var", buttons = c("fitted", 
                                                         "predictors", "other"), labels = gettextRcmdr(c("Fitted values", 
                                                                                                         "Explanatory variables", "Other (specify)")), title = gettextRcmdr("Variance Formula"), 
                 initialValue = dialog.values$initial.var)
    tkgrid(varFrame, sticky = "w")
    modelFormula(optionsFrame, hasLhs = FALSE, rhsExtras=TRUE, formulaLabel="")
    tkgrid(formulaFrame, sticky = "w")
    tkgrid(outerOperatorsFrame)
    tkgrid(getFrame(xBox), sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

DurbinWatsonTest <- function () {
    .activeModel <- ActiveModel()
	if (is.null(.activeModel)) 
		return()
	Library("lmtest")
	defaults <- list (initial.altHypothesis = "greater")
	dialog.values <- getDialog ("DurbinWatsonTest", defaults)
	initializeDialog(title = gettextRcmdr("Durbin-Waton Test"))
	tkgrid(labelRcmdr(top, text = gettextRcmdr("Test for First-Order Error Autocorrelation"), 
					fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
	onOK <- function() {
		altHypothesis <- tclvalue(altHypothesisVariable)
		putDialog ("DurbinWatsonTest", list(initial.altHypothesis = altHypothesis))
		closeDialog()
		model.formula <- as.character(formula(get(ActiveModel())))
		model.formula <- paste(model.formula[2], "~", model.formula[3])
		command <- paste("dwtest(", model.formula, ", alternative=\"", 
				altHypothesis, "\", data=", ActiveDataSet(), ")", 
				sep = "")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "dwtest", reset = "DurbinWatsonTest")
	radioButtons(name = "altHypothesis", buttons = c("greater", 
					"notequal", "less"), values = c("greater", "two.sided", 
					"less"), labels = c("rho >  0", "rho != 0", "rho <  0"), 
			title = gettextRcmdr("Alternative Hypothesis"), 
			initialValue = dialog.values$initial.altHypothesis)
	tkgrid(altHypothesisFrame, sticky = "w")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix()
}

RESETtest <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) 
        return()
    Library("lmtest")
    defaults <- list (initial.square = 1, initial.cube = 1, initial.type = "regressor")
    dialog.values <- getDialog ("RESETtest", defaults)
    initializeDialog(title = gettextRcmdr("RESET Test"))
    tkgrid(labelRcmdr(top, text = gettextRcmdr("Test for Nonlinearity"), 
                      fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    onOK <- function() {
        type <- tclvalue(typeVariable)
        square <- tclvalue(squareVariable)
        cube <- tclvalue(cubeVariable)
        putDialog ("RESETtest", list (initial.square = square, initial.cube = cube, initial.type = type))
        closeDialog()
        model.formula <- as.character(formula(get(ActiveModel())))
        model.formula <- paste(model.formula[2], "~", model.formula[3])
        if (square == "0" && cube == "0") {
            errorCondition(recall = RESETtest, message = gettextRcmdr("No powers are checked."))
            return()
        }
        powers <- if (square == "1" && cube == "1") 
            "2:3"
        else if (square == "1" && cube == "0") 
            "2"
        else if (square == "0" && cube == "1") 
            "3"
        command <- paste("resettest(", model.formula, ", power=", 
                         powers, ", type=\"", type, "\", data=", ActiveDataSet(), 
                         ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "reset", reset = "RESETtest", apply = "RESETtest")
    optionsFrame <- tkframe(top)
    squareVariable <- tclVar(dialog.values$initial.square)
    squareCheckBox <- ttkcheckbutton(optionsFrame, variable = squareVariable)
    cubeVariable <- tclVar(dialog.values$initial.cube)
    cubeCheckBox <- ttkcheckbutton(optionsFrame, variable = cubeVariable)
    typeVariable <- tclVar("regressor")
    radioButtons(optionsFrame, name = "type", buttons = c("regressor", 
                                                          "fitted", "princomp"), labels = gettextRcmdr(c("Explanatory variables", 
                                                                                                         "Fitted values", "First principal component")), title = gettextRcmdr("Type of Test"), 
                 initialValue = dialog.values$initial.type)
    tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Powers to Include"), 
                      fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("2 (squares)")), 
           squareCheckBox, sticky = "w")
    tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("3 (cubes)   ")), 
           cubeCheckBox, sticky = "w")
    tkgrid(typeFrame, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

OutlierTest <- function(){
    .activeModel <- ActiveModel()
	if (is.null(.activeModel)) return()
	if (!checkMethod("outlierTest", .activeModel)) {
		errorCondition(gettextRcmdr("There is no appropriate outlierTest method for a model of this class."))
		return()
	}
	doItAndPrint(paste("outlierTest(", .activeModel, ")", sep=""))
}

confidenceIntervals <- function () {
    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) 
        return()
    Library("MASS")
    defaults <- list (initial.level = "0.95", initial.statistic="LR")
    dialog.values <- getDialog ("confidenceIntervals", defaults)
    initializeDialog(title = gettextRcmdr("Confidence Intervals"))
    tkgrid(labelRcmdr(top, text = gettextRcmdr("Confidence Intervals for Individual Coefficients"), 
                      fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    onOK <- function() {
        level <- tclvalue(confidenceLevel)
        opts <- options(warn = -1)
        lev <- as.numeric(level)
        options(opts)
        closeDialog()
        if ((is.na(lev)) || !is.numeric(lev) || (lev < 0) || (lev > 1)) {
            Message(gettextRcmdr("Confidence level must be a number between 0 and 1."),
                    type="error")
            confidenceIntervals()
            return()
        }
        putDialog ("confidenceIntervals", list (initial.level = level,
                                                initial.statistic = if(glm) tclvalue(typeVariable) else "LR"))
        command <- if (glm) 
            paste("Confint(", .activeModel, ", level=", level, 
                  ", type=\"", tclvalue(typeVariable), "\")", sep = "")
        else paste("Confint(", .activeModel, ", level=", level, 
                   ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "Confint", reset = "confidenceIntervals", apply = "confidenceIntervals")
    confidenceFrame <- tkframe(top)
    confidenceLevel <- tclVar(dialog.values$initial.level)
    confidenceField <- ttkentry(confidenceFrame, width = "6", 
                                textvariable = confidenceLevel)
    radioButtons(top, name = "type", buttons = c("LR", "Wald"), initialValue=dialog.values$initial.statistic,
                 labels = gettextRcmdr(c("Likelihood-ratio statistic", 
                                         "Wald statistic")), title = gettextRcmdr("Test Based On"))
    tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
           confidenceField, sticky = "w")
    tkgrid(confidenceFrame, sticky = "w")
    glm <- class(get(.activeModel))[1] == "glm"
    if (glm) 
        tkgrid(typeFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

aic <- function(){
	.activeModel <- ActiveModel()
	if (is.null(.activeModel)) return()
	doItAndPrint(paste("AIC(", .activeModel, ")", sep=""))
}

bic <- function(){
	.activeModel <- ActiveModel()
	if (is.null(.activeModel)) return()
	doItAndPrint(paste("BIC(", .activeModel, ")", sep=""))
}

stepwiseRegression <- function () {
    Library("MASS")
    defaults <- list (initial.direction = "backward/forward", initial.criterion = "BIC")
    dialog.values <- getDialog ("stepwiseRegression", defaults)
    initializeDialog(title = gettextRcmdr("Stepwise Model Selection"))
    onOK <- function() {
        direction <- as.character(tclvalue(directionVariable))
        criterion <- as.character(tclvalue(criterionVariable))
        putDialog ("stepwiseRegression", list (initial.direction = tclvalue(directionVariable), 
                                               initial.criterion = tclvalue(criterionVariable)))
        closeDialog()
        doItAndPrint(paste("stepwise(", ActiveModel(), ", direction='", 
                           direction, "', criterion='", criterion, "')", sep = ""))
        tkdestroy(top)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "stepwise", reset = "stepwiseRegression", 
                 apply =  "stepwiseRegression")
    radioButtons(top, name = "direction", buttons = c("bf", "fb", 
                                                      "b", "f"), values = c("backward/forward", "forward/backward", 
                                                                            "backward", "forward"), labels = gettextRcmdr(c("backward/forward", 
                                                                                                                            "forward/backward", "backward", "forward")), title = gettextRcmdr("Direction"), 
                 initialValue = dialog.values$initial.direction)
    radioButtons(top, name = "criterion", buttons = c("bic", 
                                                      "aic"), values = c("BIC", "AIC"), labels = gettextRcmdr(c("BIC", 
                                                                                                                "AIC")), title = gettextRcmdr("Criterion"), initialValue = dialog.values$initial.criterion)
    tkgrid(directionFrame, criterionFrame, sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix()
}

subsetRegression <- function () {
    Library("leaps")
    defaults <- list (initial.criterion = "bic", initial.nbest = 1)
    dialog.values <- getDialog ("subsetRegression", defaults)
    initializeDialog(title = gettextRcmdr("Subset Model Selection"))
    onOK <- function() {
        formula <- paste(sub("^[ ]*", "", deparse(formula(get(ActiveModel())))), 
                         collapse = "")
        criterion <- as.character(tclvalue(criterionVariable))
        nbest <- as.numeric(tclvalue(nbestValue))
        putDialog ("subsetRegression", list (initial.criterion = criterion, initial.nbest = nbest))
        nvmax <- as.numeric(tclvalue(nvmaxValue))
        really.big <- if (nvmax > 50) 
            "TRUE"
        else "FALSE"
        closeDialog()
        doItAndPrint(paste("plot(regsubsets(", formula, ", data=", 
                           ActiveDataSet(), ", nbest=", nbest, ", nvmax=", nvmax, 
                           "), scale='", criterion, "')", sep = ""))
        tkdestroy(top)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "regsubsets", reset = "subsetRegression", apply = "subsetRegression")
    radioButtons(top, name = "criterion", buttons = c("bic", 
                                                      "Cp", "adjr2", "r2"), labels = gettextRcmdr(c("BIC", 
                                                                                                    "Mallows Cp", "Adjusted R-sq.", "R-squared")), title = gettextRcmdr("Criterion for Model Plot"), 
                 initialValue = dialog.values$initial.criterion)
    nvar <- ncol(model.matrix(get(ActiveModel())))
    nbestValue <- tclVar(dialog.values$initial.nbest)
    nvmaxValue <- tclVar(as.character(min(25, nvar)))
    slidersFrame <- tkframe(top)
    nbestSlider <- tkscale(slidersFrame, from = 1, to = 10, showvalue = TRUE, 
                           variable = nbestValue, resolution = 1, orient = "horizontal")
    nvmaxSlider <- tkscale(slidersFrame, from = 1, to = nvar, 
                           showvalue = TRUE, variable = nvmaxValue, resolution = 1, 
                           orient = "horizontal")
    tkgrid(tklabel(slidersFrame, text = "     "), tklabel(slidersFrame, 
                                                          text = gettextRcmdr("Number of best models\nof each size:"), 
                                                          fg = getRcmdr("title.color"), font="RcmdrTitleFont"), nbestSlider, sticky = "w")
    tkgrid(tklabel(slidersFrame, text = "     "), tklabel(slidersFrame, 
                                                          text = gettextRcmdr("Maximum size:"), fg = getRcmdr("title.color"), font="RcmdrTitleFont"), nvmaxSlider, 
           sticky = "e")
    tkgrid(criterionFrame, slidersFrame, sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix()
}

effectPlots <- function () {
  Library("effects")
  defaults <- list(initial.all.or.pick = "TRUE", initial.predictors = NULL, 
                   initial.partial.res = 0, initial.span = 50, initial.style = "stacked")
  dialog.values <- getDialog("effectPlots", defaults)
  initializeDialog(title = gettextRcmdr("Model Effect Plots"))
  predictors <- all.vars(formula(get(activeModel(), envir=.GlobalEnv))[[3]])
  predictorsFrame <- tkframe(top)
  radioButtons(predictorsFrame, name = "allEffects", buttons = c("yes", "no"), 
               values = c("TRUE", "FALSE"),  
               labels = gettextRcmdr(c("Yes", "No")), title = gettextRcmdr("Plot All High-Order Effects?"),
               initialValue = dialog.values$initial.all.or.pick)
  predictorsBox <- variableListBox(predictorsFrame, predictors, selectmode = "multiple", 
                                   title = gettextRcmdr("Predictors (pick one or more)"), 
                                   initialSelection = varPosn(dialog.values$initial.predictors, vars=predictors))
  
  partialResFrame <- tkframe(top)
  partialResVariable <- tclVar(dialog.values$initial.partial.res)
  partialResCheckBoxFrame <- tkframe(partialResFrame)
  partialResCheckBox <- ttkcheckbutton(partialResCheckBoxFrame, variable = partialResVariable)
  styleFrame <- tkframe(top)
  radioButtons(styleFrame, name="styleButtons", buttons = c("stacked", "lines"),
               labels = gettextRcmdr(c("Stacked areas", "Lines with confidence bands")),
               title=gettextRcmdr("Style of Graph"), initialValue=dialog.values$initial.style)
  sliderValue <- tclVar(dialog.values$initial.span)
  sliderFrame <- tkframe(partialResFrame)
  slider <- tkscale(sliderFrame, from = 5, to = 100, showvalue = TRUE,
                    variable = sliderValue, resolution = 5, orient = "horizontal")
  onOK <- function() {
    predictors <- getSelection(predictorsBox)
    allEffects <- as.logical(tclvalue(allEffectsVariable))
    partial.residuals <- tclvalue(partialResVariable) == "1"
    span <- as.numeric(tclvalue(sliderValue))
    style <- tclvalue(styleButtonsVariable)
    closeDialog() 
    if (allEffects){
      command <- if (class(get(activeModel(), envir=.GlobalEnv))[1] %in% c("multinom", "polr"))
        paste("plot(allEffects(", activeModel(), '), style="', style, '")', sep="")
      else paste("plot(allEffects(", activeModel(), 
                 if (partial.residuals) paste(", partial.residuals=TRUE), span=", span/100, ")", sep="")
                 else "))", sep="")
      doItAndPrint(command)
      predictors <- NULL
    }
    else {
      if (length(predictors) == 0) {
        errorCondition(recall = effectPlots, 
                       message = gettextRcmdr("You must select one or more predictors\n or plot all high-order effects."))
        return()
      }
      if (partial.residuals && (all(predictors %in% Factors()))){
        errorCondition(recall = effectPlots, 
                       message = gettextRcmdr("To plot partial residuals,\n there must be a least one numeric predictor."))
        return()
      }
      command <- if (class(get(activeModel(), envir=.GlobalEnv))[1] %in% c("multinom", "polr"))
        paste("plot(Effect(c(", paste(paste('"', predictors, '"', sep=""), collapse=", "), "), ", 
              activeModel(), '), style="', style, '")', sep="")      
      else paste("plot(Effect(c(", paste(paste('"', predictors, '"', sep=""), collapse=", "), "), ", activeModel(),
                 if (partial.residuals) paste(", partial.residuals=TRUE), span=", span/100, ")", sep="")
                 else "))", sep = "")
      doItAndPrint(command)
    }
    putDialog ("effectPlots", list(initial.all.or.pick=as.character(allEffects), initial.predictors=predictors, 
                                   initial.partial.res=as.numeric(partial.residuals),
                                   initial.span=span, initial.style=style))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "Effect", reset = "effectPlots", apply = "effectPlots")
  tkgrid(allEffectsFrame, sticky="w")
  tkgrid(getFrame(predictorsBox), sticky="w")
  tkgrid(predictorsFrame, sticky="w")
  if (class(get(activeModel(), envir=.GlobalEnv))[1] %in% c("lm", "glm")){
    tkgrid(labelRcmdr(partialResFrame, text=" "))
    tkgrid(partialResCheckBox, 
           labelRcmdr(partialResCheckBoxFrame, text=gettextRcmdr("Plot partial residuals")), 
           sticky="w")
    tkgrid(partialResCheckBoxFrame, sticky="w")
    tkgrid(slider, labelRcmdr(sliderFrame, text = gettextRcmdr("Span for smooth")),
           sticky = "swe", padx=6, pady=6)
    tkgrid(sliderFrame, sticky="w")
    tkgrid(partialResFrame, sticky="w")
  }
  else if (class(get(activeModel(), envir=.GlobalEnv))[1] %in% c("multinom", "polr")){
    tkgrid(labelRcmdr(styleFrame, text=" "))
    tkgrid(styleButtonsFrame, sticky="w")
    tkgrid(styleFrame, sticky="w")
  }
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

# effectPlots <- function(){
#   Library("effects")
#   .activeModel <- ActiveModel()
#   if (is.null(.activeModel) || !checkMethod("Effect", .activeModel)) return()
#   doItAndPrint('trellis.device(theme="col.whitebg")')
#   command <- paste("plot(allEffects(", .activeModel, "), ask=FALSE)", sep="")
#   justDoIt(command)
#   logger(command)
#   activateMenus()
#   NULL
# }
