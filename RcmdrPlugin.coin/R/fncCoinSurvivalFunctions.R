
NumericOrDate <- function(dataSet=ActiveDataSet()) {#copied from John Fox's R-RcmdrPlugin.survival, needed below 
  setdiff(Variables(), Factors())
}
startStop <- function(time){#copied from John Fox's R-RcmdrPlugin.survival, needed below
  times <- na.omit(eval(parse(text=paste(ActiveDataSet(), '[,c("', time[1], '", "', time[2],'")]', sep=""))))
  if (all(times[[time[1]]] <= times[[time[2]]])){
    return(list(start=time[1], stop=time[2], error=FALSE))
  } else if (all(times[[time[2]]] <= times[[time[1]]])){
    return(list(start=time[2], stop=time[1], error=FALSE))
  }
  else return(list(start="", stop="", error=TRUE))
}

fncCoinMaxstatSurvTest <- function(){#based on R-RcmdrPlugin.survival, survdif, and it is interlinked in this way with survival data definition
  #require(survival)
  if (!activeDataSetP()) return()
  currentModel <- FALSE
  initializeDialog(title=gettext("Maximally Selected Statistics Test for Survival data"))#, domain="R-RcmdrPlugin.survival"
  onOK <- function(){
    time <- getSelection(timeBox)
    if (length(time) == 1){
      time1 <- time
      time2 <- numeric(0)
    }
    else if (length(time) == 2){
      ss <- startStop(time)
      if (ss$error) errorCondition(recall=fncCoinMaxstatSurvTest, 
                                   message=gettext("Start and stop times must be ordered."), model=TRUE)#, domain="R-RcmdrPlugin.survival"
      time1 <- ss$start
      time2 <- ss$stop
    }
    else {
      errorCondition(recall=fncCoinMaxstatSurvTest, message=gettext("You must select one or two time variables."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    event <- getSelection(eventBox)
    if (length(event) == 0) {
      errorCondition(recall=fncCoinMaxstatSurvTest, message=gettext("You must select an event indicator."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    x <- getSelection(xBox) 
    if (length(x) == 0) {
      errorCondition(recall=fncCoinMaxstatSurvTest, message=gettext("You must select at least one explanatory variable."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    closeDialog()
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext("<all valid cases>") #, domain="R-RcmdrPlugin.survival"
        || trim.blanks(subset) == ""){
      subset <- ""
    }
    else{
      subset <- paste(", subset=", subset, sep="")
    }
    formula <- paste("Surv(", time1, ",",
                     if(length(time2) != 0) paste(time2, ",", sep=""),
                     event, ")", sep="")
    formula <- paste(formula, " ~ ", paste(x, collapse=" + "), sep="")
    command <- paste("maxstat_test(", formula, ", ", 
                     "ytrafo = function(data) trafo(data, surv_trafo = function(x) logrank_trafo(x, ties = 'HL')) ",
                     ', data=', ActiveDataSet(), subset, ")", sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="maxstat_test")
  survFrame <- tkframe(top)
  .activeDataSet <- ActiveDataSet()
  .numeric <- NumericOrDate()
  .factors <- Factors()
  time1 <- eval(parse(text=paste('attr(', .activeDataSet, ', "time1")', sep="")))
  time1 <- if (!is.null(time1)) which(time1 == .numeric) - 1 
  time2 <- eval(parse(text=paste('attr(', .activeDataSet, ', "time2")', sep="")))
  time2 <- if (!is.null(time2)) which(time2 == .numeric) - 1 
  event <- eval(parse(text=paste('attr(', .activeDataSet, ', "event")', sep="")))
  event <- if (!is.null(event)) which(event == Numeric()) - 1 
  #strata <- eval(parse(text=paste('attr(', .activeDataSet, ', "strata")', sep="")))
  #strata <- if (!is.null(strata)) which(is.element(.factors, strata)) - 1 else -1
  timeBox <- variableListBox(survFrame, NumericOrDate(), 
                             title=gettext("Time or start/end times\n(select one or two)"),
                             selectmode="multiple", initialSelection=if(is.null(time1)) NULL else c(time1, time2))
  eventBox <- variableListBox(survFrame, Numeric(), 
                              title=gettext("Event indicator\n(select one)"),
                              initialSelection=event)
  xBox <- variableListBox(survFrame, Numeric(), 
                          title=gettext("Explanatory variables\n(select one or more)"), 
                          selectmode="multiple")
  #rhoFrame <- tkframe(top)
  #rhoValue <- tclVar("0")
  #rhoSlider <- tkscale(rhoFrame, from=0, to=1, showvalue=TRUE, variable=rhoValue,
  #	resolution=0.1, orient="horizontal")
  #	modelFormula(hasLhs=FALSE)
  subsetBox()
  tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), sticky="sw")
  tkgrid(labelRcmdr(survFrame, text=""))
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(survFrame, sticky="nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(subsetFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=9, columns=1)
}

fncCoinSurvTest <- function(){#based on R-RcmdrPlugin.survival, survdif, and it is interlinked in this way with survival data definition
  #require(survival)
  if (!activeDataSetP()) return()
  currentModel <- FALSE
  initializeDialog(title=gettext("Independent Two/K Sample Test for Censored Data..."))#, domain="R-RcmdrPlugin.survival"
  onOK <- function(){
    time <- getSelection(timeBox)
    if (length(time) == 1){
      time1 <- time
      time2 <- numeric(0)
    }
    else if (length(time) == 2){
      ss <- startStop(time)
      if (ss$error) errorCondition(recall=fncCoinSurvTest, 
                                   message=gettext("Start and stop times must be ordered."), model=TRUE)#, domain="R-RcmdrPlugin.survival"
      time1 <- ss$start
      time2 <- ss$stop
    }
    else {
      errorCondition(recall=fncCoinSurvTest, message=gettext("You must select one or two time variables."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    event <- getSelection(eventBox)
    if (length(event) == 0) {
      errorCondition(recall=fncCoinSurvTest, message=gettext("You must select an event indicator."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    strata <- getSelection(strataBox) 
    if (length(strata) == 0) {
      errorCondition(recall=fncCoinSurvTest, message=gettext("You must select strata."))#, domain="R-RcmdrPlugin.survival"
      return()
    }
    block <- getSelection(blockBox) #
    if (length(block) > 0) {
      if (is.element(block, strata)) {
        errorCondition(recall=fncCoinSurvTest, message=gettextRcmdr("The group and strata variables must be different."))
        return()
      }
      block = paste(" | ", block, " ", sep="") #
      strConfintText = "" # cannot compute test for blocks !!!
    } else {
      block = "" #
    }
    Ties <- as.character(tclvalue(tiesmethodVariable)) #
    if (Ties == "default") {
      strTies = ""
    } else {
      strTies = paste(", ties.method='", Ties, "' ", sep="")
    }
    
    test <- as.character(tclvalue(testVariable))
    closeDialog()
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettext("<all valid cases>") #, domain="R-RcmdrPlugin.survival"
        || trim.blanks(subset) == ""){
      subset <- ""
    }
    else{
      subset <- paste(", subset=", subset, sep="")
    }
    
    formula <- paste("Surv(", time1, ",",
                     if(length(time2) != 0) paste(time2, ",", sep=""),
                     event, ")", sep="")
    formula <- paste(formula, " ~ ", paste(strata, collapse=" + "), sep="")
    
    if (test == "default"){
      command <- paste("surv_test(", formula, block, strTies,
                       ', data=', ActiveDataSet(), subset, ")", sep="")
    }  else command <- paste("surv_test(", formula, block, strTies,
                             ', distribution="', test, '", data=', ActiveDataSet(), subset, ")", sep="")
    
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="surv_test")
  survFrame <- tkframe(top)
  .activeDataSet <- ActiveDataSet()
  .numeric <- NumericOrDate()
  .factors <- Factors()
  time1 <- eval(parse(text=paste('attr(', .activeDataSet, ', "time1")', sep="")))
  time1 <- if (!is.null(time1)) which(time1 == .numeric) - 1 
  time2 <- eval(parse(text=paste('attr(', .activeDataSet, ', "time2")', sep="")))
  time2 <- if (!is.null(time2)) which(time2 == .numeric) - 1 
  event <- eval(parse(text=paste('attr(', .activeDataSet, ', "event")', sep="")))
  event <- if (!is.null(event)) which(event == Numeric()) - 1 
  strata <- eval(parse(text=paste('attr(', .activeDataSet, ', "strata")', sep="")))
  strata <- if (!is.null(strata)) which(is.element(.factors, strata)) - 1 else -1
  timeBox <- variableListBox(survFrame, NumericOrDate(), 
                             title=gettext("Time or start/end times\n(select one or two)"),
                             selectmode="multiple", initialSelection=if(is.null(time1)) NULL else c(time1, time2))
  eventBox <- variableListBox(survFrame, Numeric(), 
                              title=gettext("Event indicator\n(select one)"),
                              initialSelection=event)
  strataBox <- variableListBox(survFrame, Factors(), 
                               title=gettext("Strata\n(select one)"), 
                               initialSelection=strata)
  blockBox <- variableListBox(survFrame, Factors(), title="Block\n(pick none or one)") #
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="test", buttons=c("default", "exact", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Exact", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test"))       
  
  radioButtons(optionsFrame, name="tiesmethod",
               buttons=c("default","logrank", "HL", "averagescores"), #
               values=c("default","logrank", "HL", "average-scores"), initialValue="default",
               labels=gettext(c("Default","Logrank", "HL", "Average-scores")),
               title=gettext("Ties method"))
  #	modelFormula(hasLhs=FALSE)
  subsetBox()
  tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), labelRcmdr(survFrame, text="  "), getFrame(blockBox), sticky="sw")
  tkgrid(labelRcmdr(survFrame, text=""))
  tkgrid(getFrame(strataBox), sticky="nw")
  tkgrid(survFrame, sticky="nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(tiesmethodFrame, sticky="new")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(optionsFrame, sticky="nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(testFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(subsetFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=9, columns=1)
}
