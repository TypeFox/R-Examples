# last modified 2015-08-31 by J. Fox

CoxModel <- function(){
    defaults <- list(time1=NULL, time2=NULL, event=NULL, strata=NULL, cluster=NULL, 
        ties="efron", robust="default", subset=NULL, initial.tab=0, survtype=NULL)
    dialog.values <- getDialog("CoxModel", defaults)
    if (!activeDataSetP()) return()
    initializeDialog(title=gettext("Cox-Regression Model", domain="R-RcmdrPlugin.survival"),
        use.tabs=TRUE, tabs=c("dataTab", "modelTab"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "coxph"
    else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), hasLhs=TRUE)
        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
    }
    UpdateModelNumber()
    modelName <- tclVar(paste("CoxModel.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(modelTab)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        time <- getSelection(timeBox)
        if (length(time) == 1){
            time1 <- time
            time2 <- numeric(0)
        }
        else if (length(time) == 2){
            ss <- startStop(time)
            if (ss$error) errorCondition(recall=CoxModel, 
                message=gettext("Start and stop times must be ordered.", 
                    domain="R-RcmdrPlugin.survival"), model=TRUE)
            time1 <- ss$start
            time2 <- ss$stop
        }
        else {
            errorCondition(recall=CoxModel, message=gettext("You must select one or two time variables.", 
                domain="R-RcmdrPlugin.survival"), model=TRUE)
            return()
        }
        event <- getSelection(eventBox)
        survtype <- as.character(tclvalue(survtypeVariable))
        if (length(event) == 0) {
            errorCondition(recall=CoxModel, message=gettext("You must select an event indicator.", 
                domain="R-RcmdrPlugin.survival"), model=TRUE)
            return()
        }
        if (length(time) == 2 && (! survtype %in% c("default", "counting"))){
          errorCondition(recall=survregModel,
                         message=gettext("start-end times only for counting-process data.",
                                         domain="R-RcmdrPlugin.survival"))
          return()
        }
        if (length(time) == 1 && survtype == "counting"){
          errorCondition(recall=survregModel,
                         message=gettext("start-end times required for counting-process data.",
                                         domain="R-RcmdrPlugin.survival"))
          return()
        }
        strata <- getSelection(strataBox)
        cluster <- getSelection(clusterBox)
        modelValue <- trim.blanks(tclvalue(modelName))
        ties <- as.character(tclvalue(tiesVariable))
        robust <- as.character(tclvalue(robustVariable))
        subset <- tclvalue(subsetVariable)
        putDialog("CoxModel", list(
            time1=time1,
            time2=if (length(time2) == 0) NULL else time2,
            event=event, strata=strata, cluster=cluster, 
            ties=ties, robust=robust, subset=subset, initial.tab=tab, survtype=survtype
        ))
        closeDialog()
        if (!is.valid.name(modelValue)){
            errorCondition(recall=CoxModel, message=sprintf(gettext('"%s" is not a valid name.', 
                domain="R-RcmdrPlugin.survival"), modelValue), model=TRUE)
            return()
        }
        if (trim.blanks(subset) == gettext("<all valid cases>", domain="R-RcmdrPlugin.survival") 
            || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
        }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
        }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=CoxModel, message=gettext("Right-hand side of model empty.", 
                domain="R-RcmdrPlugin.survival"), model=TRUE)
            return()
        }
        if (is.element(modelValue, listCoxModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettext("Model", 
                domain="R-RcmdrPlugin.survival")))){
                UpdateModelNumber(-1)
                CoxModel()
                return()
            }
        }
        formula <- paste("Surv(", time1,
                         if (length(time2) != 0) paste(",", time2),
                         if (length(event) != 0) paste(",", event),
                         if (survtype != "default") paste(', type="', survtype, '"', sep=""),
                         ")", sep="")
        formula <- paste(formula, " ~ ", tclvalue(rhsVariable), sep="")
        if (length(strata) > 0 && length(grep("strata\\(", formula)) == 0) 
            formula <- paste(formula, " + strata(", paste(strata, collapse=","), ")", sep="")
        if (length(cluster) > 0 && length(grep("cluster\\(", formula)) == 0) 
            formula <- paste(formula, " + cluster(", cluster, ")", sep="")
        command <- paste("coxph(", formula, ', method="', ties, '"', 
            if (robust != "default") paste(", robust=", robust, sep=""),
            ", data=", ActiveDataSet(), subset, ")", sep="")
        doItAndPrint(paste(modelValue, "<-", command))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="coxph", model=TRUE, reset="CoxModel", apply="CoxModel")
    tkgrid(labelRcmdr(modelFrame, text=gettext("Enter name for model:", 
        domain="R-RcmdrPlugin.survival")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(labelRcmdr(modelTab, text=""))
    survFrame <- tkframe(dataTab)
    .activeDataSet <- ActiveDataSet()
    .numeric <- NumericOrDate()
    .factors <- Factors()
    .variables <- Variables()
    time1 <- if(!is.null(dialog.values$time1)) dialog.values$time1 else eval(parse(text=paste('attr(', .activeDataSet, ', "time1")', sep="")))
    time1 <- if (!is.null(time1)) which(time1 == .numeric) - 1 
    time2 <- if(!is.null(dialog.values$time2)) dialog.values$time2 else eval(parse(text=paste('attr(', .activeDataSet, ', "time2")', sep="")))
    time2 <- if (!is.null(time2)) which(time2 == .numeric) - 1 
    event <- if(!is.null(dialog.values$event)) dialog.values$event else eval(parse(text=paste('attr(', .activeDataSet, ', "event")', sep="")))
    event <- if (!is.null(event)) which(event == Numeric()) - 1 
    survtype <- if(!is.null(dialog.values$survtype)) dialog.values$survtype else eval(parse(text=paste('attr(', .activeDataSet, ', "survtype")', sep="")))
    if (is.null(survtype)) survtype <- "default"
    strata <- if(!is.null(dialog.values$strata)) dialog.values$strata else eval(parse(text=paste('attr(', .activeDataSet, ', "strata")', sep="")))
    strata <- if (!is.null(strata)) which(is.element(.factors, strata)) - 1 else -1
    cluster <- if(!is.null(dialog.values$cluster)) dialog.values$cluster else eval(parse(text=paste('attr(', .activeDataSet, ', "cluster")', sep="")))
    cluster <- if (!is.null(cluster)) which(cluster == if (allVarsClusters()) .variables else .factors) - 1 else -1
    timeBox <- variableListBox(survFrame, NumericOrDate(), 
        title=gettext("Time or start/end times\n(select one or two)", domain="R-RcmdrPlugin.survival"),
        selectmode="multiple", initialSelection=if(is.null(time1)) NULL else c(time1, time2))
    eventBox <- variableListBox(survFrame, Numeric(), title=gettext("Event indicator\n(select one)", 
        domain="R-RcmdrPlugin.survival"), initialSelection=event)
    strataBox <- variableListBox(survFrame, Factors(), title=gettext("Strata\n(select zero or more)", 
        domain="R-RcmdrPlugin.survival"), selectmode="multiple", initialSelection=strata)
    clusterBox <- variableListBox(survFrame, if (allVarsClusters()) Variables() else Factors(), 
        title=gettext("Clusters\n(optional)", domain="R-RcmdrPlugin.survival"), initialSelection=cluster)
    radioButtons(survFrame, name="survtype",
                 buttons=c("default", "right", "left", "counting"),
                 labels=gettext(c("Default", "Right", "Left", "Counting process")),
                 initialValue=survtype, title=gettext("Type of Censoring", domain="R-RcmdrPlugin.survival"))
    optionsFrame <- tkframe(modelTab)
    radioButtons(optionsFrame, name="ties",
        buttons=c("efron", "breslow", "exact"), initialValue=dialog.values$ties,
        labels=gettext(c("Efron", "Breslow", "Exact"), domain="R-RcmdrPlugin.survival"), 
        title=gettext("Method for Ties", domain="R-RcmdrPlugin.survival"))
    radioButtons(optionsFrame, name="robust",
        buttons=c("default", "TRUE", "FALSE"), initialValue=dialog.values$robust,
        labels=gettext(c("Default", "Yes", "No"), domain="R-RcmdrPlugin.survival"), 
        title=gettext("Robust Standard Errors", domain="R-RcmdrPlugin.survival"))
    modelFormula(modelTab, hasLhs=FALSE, rhsExtras=TRUE)
    subsetBox(dataTab, model=TRUE, subset.expression=dialog.values$subset)
    tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox),
           labelRcmdr(survFrame, text="  "), survtypeFrame, sticky="nw")
    tkgrid(labelRcmdr(survFrame, text=""))
    tkgrid(getFrame(strataBox), labelRcmdr(survFrame, text="  "), getFrame(clusterBox), sticky="nw")
    tkgrid(survFrame, sticky="w")
    tkgrid(tiesFrame, labelRcmdr(optionsFrame, text="                 "), robustFrame, sticky="new")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(labelRcmdr(modelTab, text=""))
    tkgrid(getFrame(xBox), sticky="w", columnspan=2)
    tkgrid(labelRcmdr(modelTab, text=""))
    tkgrid(labelRcmdr(outerOperatorsFrame, text="         "), operatorsFrame, sticky="w")
    tkgrid(outerOperatorsFrame, sticky="ew")
    tkgrid(formulaFrame, sticky="w")
    tkgrid(labelRcmdr(dataTab, text=""))
    tkgrid(subsetFrame, sticky="w")
    dialogSuffix(focus=rhsEntry, preventDoubleClick=TRUE, use.tabs=TRUE, grid.buttons=TRUE, 
        tabs=c("dataTab", "modelTab"), 
        tab.names=gettext("Data", "Model", domain="R-RcmdrPlugin.survival"))
}
