# Statistics Menu dialogs

# last modified 2014-07-26 by J. Fox

# Variances menu

twoVariancesFTest <- function () {
  defaults <- list(initial.groups = NULL, initial.response = NULL, initial.alternative = "two.sided", 
                   initial.confidenceLevel = ".95", initial.label=NULL, initial.tab=0)
  dialog.values <- getDialog("twoVariancesFTest", defaults)
  initializeDialog(title = gettextRcmdr("Two Variances F-Test"), use.tabs=TRUE)
  variablesFrame <- tkframe(dataTab)
  groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), 
                              title = gettextRcmdr("Groups (pick one)"), 
                              initialSelection = varPosn(dialog.values$initial.groups, "twoLevelFactor"))
  responseBox <- variableListBox(variablesFrame, Numeric(), 
                                 title = gettextRcmdr("Response Variable (pick one)"), 
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"))
  onOK <- function() {
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall = twoVariancesFTest, message = gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall = twoVariancesFTest, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    closeDialog()
    putDialog("twoVariancesFTest", list(initial.groups = group, initial.response = response, 
                                        initial.alternative = alternative, initial.confidenceLevel = level,
                                        initial.label=.groupsLabel, initial.tab=tab))
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("with(", .activeDataSet, ", tapply(", response, 
                       ", ", group, ",  var, na.rm=TRUE))", 
                       sep = ""))
    doItAndPrint(paste("var.test(", response, " ~ ", group, 
                       ", alternative='", alternative, "', conf.level=", 
                       level, ", data=", .activeDataSet, ")", sep = ""))
    tkfocus(CommanderWindow())
    #         tkdestroy(top)
  }
  OKCancelHelp(helpSubject = "var.test", reset = "twoVariancesFTest", apply = "twoVariancesFTest")
  radioButtons(optionsTab, name = "alternative", buttons = c("twosided", 
                                                             "less", "greater"), values = c("two.sided", "less", "greater"), 
               labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
                                       "Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"), 
               initialValue = dialog.values$initial.alternative,)
  confidenceFrame <- tkframe(optionsTab)
  confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
         getFrame(responseBox), sticky = "nw")
  tkgrid(variablesFrame, sticky = "w")
  groupsLabel(optionsTab, groupsBox = groupBox, initialText=dialog.values$initial.label)
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level:  "), 
                    fg = getRcmdr("title.color"), font="RcmdrTitleFont"), confidenceField, sticky = "w")
  tkgrid(alternativeFrame, sticky = "w")
  tkgrid(confidenceFrame, sticky = "w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

BartlettTest <- function () {
  defaults <- list(initial.group = NULL, initial.response = NULL)
  dialog.values <- getDialog("BartlettTest", defaults)
  initializeDialog(title = gettextRcmdr("Bartlett's Test"))
  variableFrame <- tkframe(top)
  groupBox <- variableListBox(variableFrame, Factors(), selectmode = "multiple", 
                              title = gettextRcmdr("Factors (pick one or more)"),
                              initialSelection = varPosn(dialog.values$initial.group, "factor"))
  responseBox <- variableListBox(variableFrame, Numeric(),  
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"),
                                 title = gettextRcmdr("Response Variable (pick one)"))
  onOK <- function() {
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall = BartlettTest, message = gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall = BartlettTest, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    closeDialog()
    putDialog("BartlettTest", list(initial.group = group, initial.response = response))
    .activeDataSet <- ActiveDataSet()
    if (length(group) == 1){
      doItAndPrint(paste("with(", .activeDataSet, ", tapply(",  
                         response, ", ", 
                         group, ", var, na.rm=TRUE))", sep = ""))
      doItAndPrint(paste("bartlett.test(", response, " ~ ", 
                         group, ", data=", .activeDataSet, ")", sep = ""))
    }
    else{
      command <- paste("with(", .activeDataSet, ", tapply(", response, 
                       ", list(", paste(paste(group, sep=""), collapse=", "), "), var, na.rm=TRUE))", sep="")
      doItAndPrint(command)
      doItAndPrint(paste("bartlett.test(", response, " ~ interaction(", 
                         paste(group, collapse=", "), "), data=", .activeDataSet, ")", sep = ""))
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "bartlett.test", reset = "BartlettTest", apply = "BartlettTest")
  tkgrid(getFrame(groupBox), labelRcmdr(variableFrame, text = "    "), 
         getFrame(responseBox), sticky = "nw")
  tkgrid(variableFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

LeveneTest <- function () {
  defaults <- list(initial.group = NULL, initial.response = NULL, initial.center = "median")
  dialog.values <- getDialog("LeveneTest", defaults)
  initializeDialog(title = gettextRcmdr("Levene's Test"))
  variableFrame <- tkframe(top)
  groupBox <- variableListBox(variableFrame, Factors(), selectmode = "multiple", 
                              title = gettextRcmdr("Factors (pick one or more)"),
                              initialSelection = varPosn(dialog.values$initial.group, "factor"))
  responseBox <- variableListBox(variableFrame, Numeric(), 
                                 title = gettextRcmdr("Response Variable (pick one)"),
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"))
  radioButtons(name = "center", buttons = c("median", "mean"), 
               labels = c(gettextRcmdr("median"), gettextRcmdr("mean")), 
               title = gettextRcmdr("Center"), initialValue = dialog.values$initial.center)
  onOK <- function() {
    group <- getSelection(groupBox)
    center <- as.character(tclvalue(centerVariable))
    if (length(group) == 0) {
      errorCondition(recall = LeveneTest, message = gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall = LeveneTest, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    closeDialog()
    putDialog("LeveneTest", list(initial.group = group, initial.response = response, 
                                 initial.center = center))
    .activeDataSet <- ActiveDataSet()
    if (length(group) == 1){
      doItAndPrint(paste("with(", .activeDataSet, ", tapply(", paste( 
        response, sep = ""), ", ", 
        paste(group, sep = ""), ", var, na.rm=TRUE))", sep = ""))
      doItAndPrint(paste("leveneTest(", response, " ~ ", group, ", data=",
                         .activeDataSet, ', center="', center, '")', sep=""))
    }
    else{
      command <- paste("with(", .activeDataSet, ", tapply(", response, 
                       ", list(", paste(paste(group,sep=""), collapse=", "), "), var, na.rm=TRUE))", sep="")
      doItAndPrint(command)
      doItAndPrint(paste("leveneTest(", response, " ~ ",
                         paste(group, collapse="*"), ", data=", .activeDataSet, ', center="', center, '")', sep = ""))
      
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "leveneTest", reset = "LeveneTest", apply = "LeveneTest")
  tkgrid(getFrame(groupBox), labelRcmdr(variableFrame, text = "    "), 
         getFrame(responseBox), sticky = "nw")
  tkgrid(variableFrame, sticky = "w")
  tkgrid(centerFrame, sticky = "nw")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

