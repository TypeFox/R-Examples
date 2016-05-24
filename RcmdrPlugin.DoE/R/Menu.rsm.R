## two instances of assign replaced by justDoIt
Menu.rsm <- function(response=NULL, factor.names=NULL){
    initializeDialog(title=gettextRcmdr("Response surface model"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
        "rsm" %in% class(get(.activeModel, envir=.GlobalEnv))[1]
        else FALSE
    
    ## generally do not use current model!
    currentModel <- FALSE
    
    oldwarn <- options("warn")$warn
    options(warn=0)
    if (is.na(as.numeric(getRcmdr("degree")))) degree <- NULL 
        else degree <- as.numeric(getRcmdr("degree"))
    options(warn=oldwarn)

    ## working with current models does not work because of 
    ## switching between coded and uncoded designs
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
        if (!identical(currentFields$data, ActiveDataSet())) currentModel <- FALSE
        }
    if (!currentModel) {
        if (!as.logical(getRcmdr("coded"))) 
            hilf <- rsmformula(get(ActiveDataSet(), envir=.GlobalEnv), degree=degree, 
            response=response, factor.names=factor.names, coded=FALSE)
        else
        hilf <- rsmformula(get(ActiveDataSet(), envir=.GlobalEnv), degree=degree, 
            response=response, factor.names=factor.names)
        currentFields <- list(lhs=as.character(hilf[2]), rhs=as.character(hilf[3]), 
           data=ActiveDataSet(), subset="")
        currentModel <- TRUE
    }
    UpdateModelNumber()
    modelName <- tclVar(paste("rsmModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)


    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()
        if (!is.valid.name(modelValue)){
            errorCondition(recall=Menu.rsm, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=Menu.rsm, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=Menu.rsm, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
            return()
            }
#        if (is.element(modelValue, listRSMs())) {
#            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
#                UpdateModelNumber(-1)
#                rsmModel()
#                return()
#                }
#            }
        formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
        if (as.logical(getRcmdr("coded"))) {
            hilfnam <- paste(getRcmdr(".activeDataSet"), "coded",sep=".")
            if (exists(hilfnam, , envir=.GlobalEnv)){
              if ("no" == tclvalue(checkReplace(hilfnam, gettextRcmdr("Object")))){
                 errorCondition(window=top,recall=Menu.rsm, 
                     message=gettextRcmdr("Execution stopped by user:\ncoded design would have overwritten existing design"))
                     return()}
            }
            ## replace assign with justDoIt; assign(hilfnam, code.design(get(ActiveDataSet())), envir=.GlobalEnv)
            cmd <- paste(hilfnam, " <- code.design(", getRcmdr(".activeDataSet"), ")")
            justDoIt(cmd)
            logger(cmd)
            putRcmdr(".activeDataSet", hilfnam)
            ## refresh active data set
            ActiveDataSet(hilfnam)
            activeDataSet(hilfnam)
        }
        command <- paste("rsm(", formula,
            ", data=", ActiveDataSet(), ")", sep="")

        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=top,recall=Menu.rsm, message=gettextRcmdr(hilf))
             return()
            }
        
        logger(paste(modelValue, " <- ", command, sep=""))
        ## replacel assign with justDoIt; assign(modelValue, hilf, envir=.GlobalEnv)
        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(modelValue, hilf, envir=.GlobalEnv)
        justDoIt(paste(modelValue, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.rsm", model=TRUE)
    tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    rsmFormula()          ## this grids the full model functionality
                            ## it is a macro
                            
    ## subsetBox(model=TRUE)   ## omit subset box? or not ?
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
    ## tkgrid(subsetFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }
