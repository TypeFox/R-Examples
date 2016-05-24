llinear0 <-function () 
{
Library("MASS")
Library("nnet")
Library("car")
Library("effects")


    initializeDialog(title = gettextRcmdr(paste("Mod","\U00E8","le lin", "\U00E9", "aire", 
            sep = "")))
    .activeModel <- ActiveModel()
    .activeDataSet <- ActiveDataSet()
    currentModel <- if (!is.null(.activeModel)) 
        {(class(get(.activeModel, envir = .GlobalEnv))[1]=="polr") |
(class(get(.activeModel, envir = .GlobalEnv))[1]=="lm") |
(class(get(.activeModel, envir = .GlobalEnv))[1]=="glm") |
(class(get(.activeModel, envir = .GlobalEnv))[1]=="multinom")
}
    else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir = .GlobalEnv))
   #     if (currentFields$data != .activeDataSet) 
    #        currentModel <- FALSE
    }
    UpdateModelNumber()
    modelName <- tclVar(paste("LModel.", getRcmdr("modelNumber"), 
        sep = ""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
    onOK <- function() {
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()
        if (!is.valid.name(modelValue)) {
            errorCondition(recall = llinear0, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
                modelValue), model = TRUE)
            return()
        }
        subset <- tclvalue(subsetVariable)
        if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || 
            trim.blanks(subset) == "") {
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
        }
        else {
            subset <- paste(", subset=", subset, sep = "")
            putRcmdr("modelWithSubset", TRUE)
        }
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall = llinear0, message = gettextRcmdr("Left-hand side of model empty."), 
                model = TRUE)
            return()
        }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall = llinear0, message = gettextRcmdr("Right-hand side of model empty."), 
                model = TRUE)
            return()
        }
        if (is.binary(eval(parse(text = tclvalue(lhsVariable)), 
            envir = get(.activeDataSet, envir = .GlobalEnv)))) {
            if (is.element(modelValue, listGeneralizedLinearModels())) {
                if ("no" == tclvalue(checkReplace(modelValue, 
                  type = gettextRcmdr("Model")))) {
                  UpdateModelNumber(-1)
                  generalizedLinearModel()
                  return()
                }
            }
            formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), 
                sep = " ~ ")
            command <- paste("glm(", formula, ", data=", .activeDataSet, 
                subset, ",family=binomial)", sep = "")
        }
        else {
            if (is.ordered(eval(parse(text = tclvalue(lhsVariable)), 
                envir = get(.activeDataSet, envir = .GlobalEnv)))) {
                if (is.element(modelValue, listProportionalOddsModels())) {
                  if ("no" == tclvalue(checkReplace(modelValue, 
                    type = gettextRcmdr("Model")))) {
                    UpdateModelNumber(-1)
                    proportionalOddsModel()
                    return()
                  }
                }
                formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), 
                  sep = " ~ ")
                command <- paste("polr(", formula, ", method=\"", 
                  "logistic", "\", data=", .activeDataSet, subset, 
                  ", Hess=TRUE)", sep = "")
            }
            else {
                if (is.factor(eval(parse(text = tclvalue(lhsVariable)), 
                  envir = get(.activeDataSet, envir = .GlobalEnv)))) {
                  if (is.element(modelValue, listMultinomialLogitModels())) {
                    if ("no" == tclvalue(checkReplace(modelValue, 
                      type = gettextRcmdr("Model")))) {
                      UpdateModelNumber(-1)
                      multinomialLogitModel()
                      return()
                    }
                  }
                  formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), 
                    sep = " ~ ")
                  command <- paste("multinom(", formula, ",data=na.omit(", 
                    .activeDataSet,")", subset, ")", sep = "")
                }
                else {
                  if (is.numeric(eval(parse(text = tclvalue(lhsVariable)), 
                    envir = get(.activeDataSet, envir = .GlobalEnv)))) {
                    if (is.element(modelValue, listLinearModels())) {
                      if ("no" == tclvalue(checkReplace(modelValue, 
                        type = gettextRcmdr("Model")))) {
                        UpdateModelNumber(-1)
                        linearModel()
                        return()
                      }
                    }
                    formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), 
                      sep = " ~ ")
                    command <- paste("lm(", formula, ",data=", 
                      .activeDataSet, subset, ")", sep = "")
                  }
                }
            }
        }

doItAndPrint(paste(modelValue, " <- ", command, sep = ""))
doItAndPrint(paste("sortAnova(", modelValue, ")", sep = ""))

        activeModel(modelValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "llinear0", model = TRUE,reset="llinear0")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model:")), 
        model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")
    modelFormula()
    subsetBox(model = TRUE)
    tkgrid(getFrame(xBox), sticky = "w")
    tkgrid(outerOperatorsFrame, sticky = "w")
    tkgrid(formulaFrame, sticky = "w")
    tkgrid(subsetFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 7, columns = 1, focus = lhsEntry, preventDoubleClick = TRUE)
}

