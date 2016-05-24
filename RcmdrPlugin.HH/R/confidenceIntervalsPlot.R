"confidenceIntervalsPlot" <-
function(){
  initializeDialog(title=gettextRcmdr("Confidence Intervals in Simple Linear Regression"))
    variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric,
        title=gettextRcmdr("Explanatory variables (pick one)"))
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Response variable (pick one)"))
    UpdateModelNumber()
    modelName <- tclVar(paste("RegModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=modelName)
    subsetBox()
    onOK <- function(){
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            UpdateModelNumber(-1)
            errorCondition(recall=confidenceIntervalsPlot, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=confidenceIntervalsPlot, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }        
        if (is.element(y, x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=confidenceIntervalsPlot, message=gettextRcmdr("Response and explanatory variables must be different."))
            return()
            }
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
        modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            UpdateModelNumber(-1)
            errorCondition(recall=confidenceIntervalsPlot, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                confidenceIntervalsPlot()
                return()
                }
            }
        command <- paste("lm(", y, "~", paste(x, collapse="+"),
            ", data=", ActiveDataSet(), subset, ")", sep="")
        justDoIt(paste(modelValue, " <- ", command, sep=""))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        doItAndPrint(paste("ci.plot(", modelValue, ")", sep=""))
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ci.plot", model=TRUE)
    tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(getFrame(yBox), tklabel(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")    
    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, stick="w")
    tkgrid.configure(helpButton, sticky="e")
    dialogSuffix(rows=4, columns=1)
    }

