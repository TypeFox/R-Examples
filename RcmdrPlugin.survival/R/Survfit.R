# last modified 2015-08-31 by J. Fox

Survfit <- function(){
    defaults <- list(time1=NULL, time2=NULL, event=NULL, strata=NULL, type="kaplan-meier", error="greenwood", 
        survtype=NULL, detail="default", conftype="log", conf.int="default", lev=".95", markTime="1", 
        quantiles=".25, .5, .75", subset=NULL, initial.tab=0)
    dialog.values <- getDialog("Survfit", defaults)
    if (!activeDataSetP()) return()
    currentModel <- FALSE
    initializeDialog(title=gettext("Survival Function", domain="R-RcmdrPlugin.survival"),
        use.tabs=TRUE)
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        time <- getSelection(timeBox)
        if (length(time) == 1){
            time1 <- time
            time2 <- numeric(0)
        }
        else if (length(time) == 2){
            ss <- startStop(time)
            if (ss$error) errorCondition(recall=Survfit, 
                message=gettext("Start and stop times must be ordered.", 
                    domain="R-RcmdrPlugin.survival"), model=TRUE)
            time1 <- ss$start
            time2 <- ss$stop
        }
        else {
            errorCondition(recall=Survfit, 
                message=gettext("You must select one or two time variables.", 
                    domain="R-RcmdrPlugin.survival"))
            return()
        }
        event <- getSelection(eventBox)
        strata <- getSelection(strataBox)
        type <- as.character(tclvalue(typeVariable))
        error <- as.character(tclvalue(errorVariable))
        survtype <- as.character(tclvalue(survtypeVariable))
        detail <- as.character(tclvalue(detailVariable))
        conftype <- as.character(tclvalue(conftypeVariable))
        conf.int <- as.character(tclvalue(plotconfVariable))
        lev <- as.numeric(tclvalue(confidenceLevel))
        markTime <- if (tclvalue(markTimeValue)  == 1) "TRUE" else "FALSE"
        quantiles <- tclvalue(quantiles)
        subset <- tclvalue(subsetVariable)
        putDialog("Survfit", list(time1=time1,
            time2=if (length(time2) == 0) NULL else time2,
            event=event, strata=strata, type=type, detail=detail, error=error, survtype=survtype, conftype=conftype, 
            conf.int=conf.int, lev=lev, markTime=if (markTime) "1" else "0", quantiles=quantiles, subset=subset,
            initial.tab=tab))
        if (survtype == "interval" && length(event) == 0){
            errorCondition(recall=Survfit, 
                message=gettext("You must select an event indicator if censoring is 'interval'.", 
                    domain="R-RcmdrPlugin.survival"))
            return()
        }
        if (survtype == "interval2" && length(event) != 0){
            errorCondition(recall=Survfit, 
                message=gettext("You should not select an event indicator if censoring is 'interval2'.", 
                    domain="R-RcmdrPlugin.survival"))
            return()
        }
        if (length(time) == 2 && (! survtype %in% c("counting", "interval", "interval2"))){
            errorCondition(recall=Survfit,
                message=gettext("start-end times only for counting-process or interval censoring.",
                    domain="R-RcmdrPlugin.survival"))
            return()
        }
        if (length(time) == 1 && survtype %in% c("counting", "interval", "interval2")){
            errorCondition(recall=Survfit,
                message=gettext("start-end times required for counting-process or interval censoring.",
                    domain="R-RcmdrPlugin.survival"))
            return()
        }
        quants <- paste("c(", gsub(",+", ",", gsub(" ", ",", quantiles)), ")", sep="")
        closeDialog()
        if ((is.na(lev)) || (lev < 0) || (lev > 1)) {
            errorCondition(recall=Survfit, message=gettext("Confidence level must be a number between 0 and 1.", 
                domain="R-RcmdrPlugin.survival"))
            return()
        }
        if (trim.blanks(subset) == gettext("<all valid cases>", domain="R-RcmdrPlugin.survival") 
            || trim.blanks(subset) == ""){
            subset <- ""
        }
        else{
            subset <- paste(", subset=", subset, sep="")
        }
        formula <- paste("Surv(", time1,
            if (length(time2) != 0) paste(",", time2),
            if (length(event) != 0) paste(",", event),
            if (survtype != "default") paste(', type="', survtype, '"', sep=""),
            ")", sep="")
        formula <- if (length(strata) > 0) paste(formula, " ~ ", paste(strata, collapse=" + "), sep="")
        else paste(formula, "~ 1")
        command <- paste("survfit(", formula, ', conf.type="', conftype, 
            '", conf.int=', lev, ', type="', type, '", error="', error,
            '", data=', ActiveDataSet(), subset, ")", sep="")
        doItAndPrint(paste(".Survfit <-", command))
        if (detail == "detailed") doItAndPrint("summary(.Survfit)") else doItAndPrint(".Survfit")
        conf.int <- if (conf.int == "default") "" else paste(", conf.int=", conf.int, sep="") 
        if (length(strata) == 0) doItAndPrint(paste("plot(.Survfit", conf.int, ", mark.time=", 
            markTime,  ")", sep=""))
        else{
            allstrata <- eval(parse(text=paste("with(", ActiveDataSet(), 
                ", interaction(", paste(strata, collapse=","), "))")))
            levels <- levels(allstrata)
            nlevels <- length(levels)
            doItAndPrint(paste("plot(.Survfit, col=1:", nlevels,", lty=1:", nlevels, 
                conf.int, ", mark.time=", markTime, ')', sep=""))
            doItAndPrint(paste('legend("bottomleft", legend=c(', paste(paste('"', levels, '"', sep=""), collapse=","), 
                '), title="', as.character(formula(.Survfit)[3]),
                '", col=1:', nlevels,', lty=1:', nlevels, ', bty="n")', sep=""))
        }
        doItAndPrint(paste("quantile(.Survfit, quantiles=", quants, ")", sep=""))
        logger("remove(.Survfit)")
        remove(.Survfit, envir=.GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="survfit", reset="Survfit", apply="Survfit")
    survFrame <- ttkframe(dataTab)
    optionsFrame <- ttkframe(optionsTab)
    .activeDataSet <- ActiveDataSet()
    .numeric <- NumericOrDate()
    .factors <- Factors()
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
    timeBox <- variableListBox(survFrame, NumericOrDate(), 
        title=gettext("Time or start/end times\n(select one or two)", domain="R-RcmdrPlugin.survival"),
        selectmode="multiple", initialSelection=if(is.null(time1)) NULL else c(time1, time2))
    eventBox <- variableListBox(survFrame, Numeric(), 
        title=gettext("Event indicator\n(select one or none)", domain="R-RcmdrPlugin.survival"),
        initialSelection=event)
    strataBox <- variableListBox(survFrame, Factors(), 
        title=gettext("Strata\n(select zero or more)", domain="R-RcmdrPlugin.survival"), 
        selectmode="multiple", initialSelection=strata)
    radioButtons(survFrame, name="survtype",
        buttons=c("default", "right", "left", "interval", "interval2", "counting"),
        labels=gettext(c("Default", "Right", "Left", "Interval", "Interval type 2", "Counting process"), domain="R-RcmdrPlugin.survival"),
        initialValue=survtype, title=gettext("Type of Censoring", domain="R-RcmdrPlugin.survival"))
    radioButtons(optionsTab, name="detail", 
        buttons=c("default", "detailed"),
        labels=gettext(c("Default", "Detailed")),
        initialValue=dialog.values$detail, title=gettext("Summary", domain="R-RcmdrPlugin.survival"))
    confidenceFrame <- ttkframe(optionsTab)
    radioButtons(confidenceFrame, name="conftype",
        buttons=c("log", "loglog", "plain", "none"), 
        values=c("log", "log-log", "plain", "none"), initialValue=dialog.values$conftype,
        labels=gettext(c("Log", "Log-log","Plain", "None", domain="R-RcmdrPlugin.survival")), 
        title=gettext("Confidence Intervals", domain="R-RcmdrPlugin.survival"))
    confidenceLevel <- tclVar(dialog.values$lev)
    confidenceFieldFrame <- ttkframe(confidenceFrame)
    confidenceField <- ttkentry(confidenceFieldFrame, width="6", textvariable=confidenceLevel)
    radioButtons(confidenceFrame, name="plotconf",
        buttons=c("default", "yes", "no"), 
        values=c("default", "TRUE", "FALSE"), initialValue=dialog.values$conf.int,
        labels=gettext(c("Default", "Yes", "No"), domain="R-RcmdrPlugin.survival"), 
        title=gettext("Plot confidence Intervals", domain="R-RcmdrPlugin.survival"))
    options2Frame <- ttkframe(optionsTab)
    radioButtons(options2Frame, name="type",
        buttons=c("kaplanmeier","flemingharrington", "fh2"), 
        values=c("kaplan-meier","fleming-harrington", "fh2"), initialValue=dialog.values$type,
        labels=gettext(c("Kaplan-Meier", "Fleming-Harrington", "Fleming-Harrington 2"), 
            domain="R-RcmdrPlugin.survival"),
        title=gettext("Method", domain="R-RcmdrPlugin.survival"))
    radioButtons(options2Frame, name="error",
        buttons=c("greenwood", "tsiatis"), initialValue=dialog.values$error,
        labels=gettext(c("Greenwood", "Tsiatis"), domain="R-RcmdrPlugin.survival"), 
        title=gettext("Variance Method", domain="R-RcmdrPlugin.survival"))
    quantilesFrame <- ttkframe(options2Frame)
    quantilesVariable <- tclVar("1")
    quantiles <- tclVar(dialog.values$quantiles)
    quantilesEntry <- ttkentry(quantilesFrame, width="20", textvariable=quantiles)
    markTimeFrame <- ttkframe(optionsTab)
    markTimeBox <- ttkcheckbutton(markTimeFrame)
    markTimeValue <- tclVar(dialog.values$markTime)
    tkconfigure(markTimeBox, variable=markTimeValue)	
    subsetBox(window=dataTab, subset.expression=dialog.values$subset)
    tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), sticky="sw")
    tkgrid(labelRcmdr(survFrame, text=""))
    tkgrid(getFrame(strataBox), labelRcmdr(survFrame, text="  "), survtypeFrame, labelRcmdr(survFrame, text="  "), sticky="nw")
    tkgrid(survFrame, sticky="nw")
    tkgrid(labelRcmdr(confidenceFieldFrame, text=gettext("Confidence level", 
        domain="R-RcmdrPlugin.survival"), foreground="blue"), sticky="nw")
    tkgrid(confidenceField, sticky="nw")
    tkgrid(conftypeFrame, labelRcmdr(confidenceFrame, text="          "), plotconfFrame, 
        labelRcmdr(confidenceFrame, text="      "), confidenceFieldFrame, sticky="nw")
    tkgrid(labelRcmdr(optionsTab, text=""))
    tkgrid(confidenceFrame, sticky="nw")
    tkgrid(labelRcmdr(quantilesFrame, text=gettext("Quantiles to estimate", 
        domain="R-RcmdrPlugin.survival"), foreground="blue"), sticky="nw")
    tkgrid(quantilesEntry)
    tkgrid(typeFrame, labelRcmdr(options2Frame, text="  "), errorFrame, 
        labelRcmdr(options2Frame, text="         "), quantilesFrame, sticky="new")
    tkgrid(labelRcmdr(optionsTab, text=""))
    tkgrid(optionsFrame, sticky="nw")
    tkgrid(options2Frame, sticky="nw")
    tkgrid(labelRcmdr(optionsTab, text=""))
    tkgrid(labelRcmdr(markTimeFrame, text=gettext("Mark censoring times",
        domain="R-RcmdrPlugin.survival"), fg="blue"), markTimeBox, sticky="nw")
    tkgrid(markTimeFrame, sticky="w")
    tkgrid(detailFrame, sticky="w")
    tkgrid(labelRcmdr(dataTab, text=""))
    tkgrid(subsetFrame, sticky="w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}
