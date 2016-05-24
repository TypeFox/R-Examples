"bestSubsetsRegressionModel.HH" <-
function(){
    initializeDialog(title=gettextRcmdr("Best Subsets Regression"))
    variablesFrame <- tkframe(top)
    .numeric <- Numeric()
    xBox <- variableListBox(variablesFrame, .numeric, selectmode="multiple",
        title=gettextRcmdr("Explanatory variables (pick one or more)"), listHeight=7)
    yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Response variable (pick one)"), listHeight=7)
    UpdateModelNumber()
    modelName <- tclVar(paste("RegModel.", getRcmdr("modelNumber"), sep=""))
    subsetsName <- tclVar(paste("Subsets.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- tkentry(modelFrame, width="20", textvariable=subsetsName)
    subsetBox()
    nbestName <- tclVar("2")
    nbestFrame <- tkframe(top)
    nbest <- tkentry(nbestFrame, width="3", textvariable=nbestName)

    radioButtons(name="statistic",
                 buttons=c("rsq", "rss", "adjr2", "cp", "bic", "stderr"),
                 values=c("rsq", "rss", "adjr2", "cp", "bic", "stderr"),
                 initialValue="adjr2",
                 labels=gettextRcmdr(c("R Square", "Residual Sum of Squares",
                   "Adjusted R^2", "Cp", "BIC", "Standard Error")),
                 title=gettextRcmdr("Statistic to plot"))

    onOK <- function() {
        x <- getSelection(xBox)
        y <- getSelection(yBox)
        closeDialog()
        if (0 == length(y)) {
            UpdateModelNumber(-1)
            errorCondition(recall=bestSubsetsRegressionModel.HH, message=gettextRcmdr("You must select a response variable."))
            return()
            }
        if (0 == length(x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=bestSubsetsRegressionModel.HH, message=gettextRcmdr("No explanatory variables selected."))
            return()
            }
        if (is.element(y, x)) {
            UpdateModelNumber(-1)
            errorCondition(recall=bestSubsetsRegressionModel.HH, message=gettextRcmdr("Response and explanatory variables must be different."))
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
        subsetsValue <- trim.blanks(tclvalue(subsetsName))
        if (!is.valid.name(subsetsValue)){
            UpdateModelNumber(-1)
            errorCondition(recall=bestSubsetsRegressionModel.HH, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), subsetsValue))
            return()
            }
        if (is.element(subsetsValue, listLinearModels())) {
            if ("no" == tclvalue(checkReplace(subsetsValue, type=gettextRcmdr("Model")))) {
                UpdateModelNumber(-1)
                bestSubsetsRegressionModel.HH()
                return()
                }
            }

        nbestValue <- as.integer(tclvalue(nbestName))
        if (!((is.integer(nbestValue)) && nbestValue>0 && nbestValue<=length(x))) {
            UpdateModelNumber(-1)
            errorCondition(recall=bestSubsetsRegressionModel.HH, message=sprintf(gettextRcmdr('"%s" is not a valid nbest value.'), nbestValue))
            return()
            }
        statisticValue <- tclvalue(statisticVariable)

        command <- paste("leaps::regsubsets(", y, "~", paste(x, collapse="+"),
                         ", data=", ActiveDataSet(), subset,
                         ", nbest=", nbestValue,
                         ")", sep="")
        doItAndPrint(paste(subsetsValue, " <- ", command, sep=""))

        command <- paste("summaryHH(", subsetsValue, ")", sep="")
        summaryValue <- paste(subsetsValue, "Summary", sep=".")
        justDoIt(paste(summaryValue, " <- ", command, sep=""))
        ## summaries <-
        doItAndPrint(summaryValue)

        command <- paste("plot(", summaryValue, ", statistic='", statisticValue, "', legend=FALSE)", sep="")
        logger(command)
        justDoIt(command)

##        bringToTop()
        if (version$os == "mingw32") justDoIt("bringToTop()")
        .nmax <- attr(get(summaryValue, envir=.GlobalEnv), "n.max.adjr2")
        ## .nmax <- attr(summaries, "n.max.adjr2")
        modelValue <- paste(trim.blanks(tclvalue(modelName)), .nmax, sep=".")

        command <- paste("lm.regsubsets(", subsetsValue, ", ", .nmax, ")", sep="")
        justDoIt(paste(modelValue, " <- ", command, "  ## subset ", .nmax, " has largest adjr2", sep=""))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))

        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="regsubsets", model=TRUE)
    tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(getFrame(yBox), tklabel(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")

    tkgrid(tklabel(nbestFrame, text=gettextRcmdr("Number of subsets of each size to record:")), nbest, sticky="w")
    tkgrid(nbestFrame, sticky="w")

    tkgrid(statisticFrame, sticky="w")

    tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, stick="w")
    tkgrid.configure(helpButton, sticky="e")
    dialogSuffix(rows=5, columns=1)
    }

