## one instance of assign replaced by justDoIt

Menu.linearModelDesign <- function(response=NULL){
    initializeDialog(title=gettextRcmdr("Linear model for experimental designs"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
        "lm" %in% class(get(.activeModel, envir=.GlobalEnv))[1]
        else FALSE
    currentModel <- FALSE
    oldwarn <- options("warn")$warn
    options(warn=-1)  ## suppress frequent harmless warning
    if (is.na(as.numeric(getRcmdr("degree")))){ 
             degree <- 2 
             if (!any(grepl("^",formula(get(ActiveDataSet())), fixed=TRUE))) degree <- 1
       }
        else degree <- as.numeric(getRcmdr("degree"))
    options(warn=oldwarn)

    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
        }
    if (!currentModel) {
        hilf <- formula(get(ActiveDataSet(), envir=.GlobalEnv), degree=degree, 
            response=response)
        currentFields <- list(lhs=as.character(hilf[2]), rhs=as.character(hilf[3]), 
           data=ActiveDataSet(), subset="")
        currentModel <- TRUE
    }
    UpdateModelNumber()
    modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)


    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()
        if (!is.valid.name(modelValue)){
            errorCondition(recall=Menu.linearModelDesign, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=Menu.linearModelDesign, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=Menu.linearModelDesign, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
            return()
            }
        if (is.element(modelValue, listLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                linearModel()
                return()
                }
            }
        formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
        command <- paste("lm(", formula,
            ", data=", ActiveDataSet(), ")", sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.linearModelDesign, message=gettextRcmdr(hilf))
             return()
            }
        cmd <- paste(modelValue, " <- ", command, sep="")
        logger(cmd)
        ## replace assign by justDoIt; assign(modelValue, hilf, envir=.GlobalEnv)
        putRcmdr("hilf", hilf)
        justDoIt(paste(modelValue, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        putRcmdr("modelWithSubset", FALSE)  ## for avoiding problems with adding observation statistics
          closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.linearModelDesign", model=TRUE)
    tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    modelFormula()          ## this grids the full model functionality
                            ## it is a macro
                            
    ## subsetBox(model=TRUE)   ## omit subset box? or not ?
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
    ## tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }
