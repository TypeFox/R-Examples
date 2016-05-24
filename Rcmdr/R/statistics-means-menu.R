# Statistics Menu dialogs

# last modified 2015-12-16 by J. Fox

# Means menu

independentSamplesTTest <- function () {
    defaults <- list(initial.group = NULL, initial.response = NULL, initial.alternative = "two.sided", 
                     initial.confidenceLevel = ".95", initial.variances = "FALSE", initial.label=NULL,
                     initial.tab=0)
    dialog.values <- getDialog("independentSamplesTTest", defaults)
    initializeDialog(title = gettextRcmdr("Independent Samples t-Test"), use.tabs=TRUE)
    variablesFrame <- tkframe(dataTab)
    groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), 
                                title = gettextRcmdr("Groups (pick one)"), 
                                initialSelection = varPosn(dialog.values$initial.group, "twoLevelFactor"))
    responseBox <- variableListBox(variablesFrame, Numeric(), 
                                   title = gettextRcmdr("Response Variable (pick one)"),
                                   initialSelection = varPosn(dialog.values$initial.response, "numeric"))
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        group <- getSelection(groupBox)
        if (length(group) == 0) {
            errorCondition(recall = independentSamplesTTest, 
                           message = gettextRcmdr("You must select a groups variable."))
            return()
        }
        response <- getSelection(responseBox)
        if (length(response) == 0) {
            errorCondition(recall = independentSamplesTTest, 
                           message = gettextRcmdr("You must select a response variable."))
            return()
        }
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceLevel)
        variances <- as.character(tclvalue(variancesVariable))
        putDialog ("independentSamplesTTest", list (initial.group = group, initial.response = response, initial.alternative = alternative, 
                                                    initial.confidenceLevel = level, initial.variances = variances, 
                                                    initial.label=.groupsLabel, initial.tab=tab))        
        closeDialog()
        doItAndPrint(paste("t.test(", response, "~", group, ", alternative='", 
                           alternative, "', conf.level=", level, ", var.equal=", 
                           variances, ", data=", ActiveDataSet(), ")", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "t.test", reset = "independentSamplesTTest", apply = "independentSamplesTTest")
    optionsFrame <- tkframe(optionsTab)
    radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
                                                                 "less", "greater"), values = c("two.sided", "less", "greater"), 
                 labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
                                         "Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"),
                 initialValue = dialog.values$initial.alternative)
    confidenceFrame <- tkframe(optionsFrame)
    confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
    confidenceField <- ttkentry(confidenceFrame, width = "6", 
                                textvariable = confidenceLevel)
    radioButtons(optionsFrame, name = "variances", buttons = c("yes", 
                                                               "no"), values = c("TRUE", "FALSE"),  
                 labels = gettextRcmdr(c("Yes", "No")), title = gettextRcmdr("Assume equal variances?"),
                 initialValue = dialog.values$initial.variances)
    tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
           getFrame(responseBox), sticky = "nw")
    tkgrid(variablesFrame, sticky = "nw")
    tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
                      fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
    tkgrid(confidenceField, sticky = "w")
    groupsLabel(optionsTab, groupsBox = groupBox, initialText=dialog.values$initial.label)
    tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text = "    "), 
           confidenceFrame, labelRcmdr(optionsFrame, text = "    "), 
           variancesFrame, sticky = "nw")
    tkgrid(optionsFrame, sticky = "nw")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

pairedTTest <- function () {
  defaults <- list(initial.x = NULL, initial.y = NULL, initial.alternative = "two.sided", 
                   initial.confidenceLevel = ".95", initial.tab=0)
  dialog.values <- getDialog("pairedTTest", defaults)
  initializeDialog(title = gettextRcmdr("Paired t-Test"), use.tabs=TRUE)
  .numeric <- Numeric()
  dataFrame <- tkframe(dataTab)
  xBox <- variableListBox(dataFrame, .numeric, title = gettextRcmdr("First variable (pick one)"),
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"))
  yBox <- variableListBox(dataFrame, .numeric, title = gettextRcmdr("Second variable (pick one)"),
                          initialSelection = varPosn(dialog.values$initial.y, "numeric"))
  onOK <- function() {
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    if (length(x) == 0 | length(y) == 0) {
      errorCondition(recall = pairedTTest, message = gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y) {
      errorCondition(recall = pairedTTest, message = gettextRcmdr("Variables must be different."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    putDialog ("pairedTTest", list (initial.x = x, initial.y = y, initial.alternative = alternative, 
                                    initial.confidenceLevel = level, initial.tab=tab))
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("with(", ActiveDataSet (), ", (t.test(", x, 
                       ", ", y, ", alternative='", 
                       alternative, "', conf.level=", level, ", paired=TRUE)))", 
                       sep = ""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "t.test", reset = "pairedTTest", apply = "pairedTTest")
  optionsFrame <- tkframe(optionsTab)
  radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
                                                               "less", "greater"), values = c("two.sided", "less", "greater"), 
               labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
                                       "Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"), 
               initialValue = dialog.values$initial.alternative)
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)
  tkgrid(getFrame(xBox), labelRcmdr(dataFrame, text="  "), getFrame(yBox), sticky = "nw")
  tkgrid(dataFrame, sticky="w")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level"), 
                    fg = getRcmdr("title.color"), font="RcmdrTitleFont"))
  tkgrid(confidenceField, sticky = "w")
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="  "), confidenceFrame, sticky = "nw")
  tkgrid(optionsFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

singleSampleTTest <- function () {
  defaults <- list (initial.x = NULL, initial.alternative = "two.sided", initial.level = ".95", 
                    initial.mu = "0.0")
  dialog.values <- getDialog ("singleSampleTTest", defaults)  
  initializeDialog(title = gettextRcmdr("Single-Sample t-Test"))
  xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"),
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"))
  onOK <- function() {
    x <- getSelection(xBox)
    if (length(x) == 0) {
      errorCondition(recall = singleSampleTTest, message = gettextRcmdr("You must select a variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    mu <- tclvalue(muVariable)
    putDialog ("singleSampleTTest", list (initial.x = x, initial.alternative = alternative, 
                                          initial.level = level, initial.mu = mu))
    closeDialog()
    doItAndPrint(paste("with(", ActiveDataSet (), ", (t.test(", x, 
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", 
                       level, ")))", sep = ""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "t.test", reset = "singleSampleTTest", apply = "singleSampleTTest")
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
                                                               "less", "greater"), values = c("two.sided", "less", "greater"), 
               labels = gettextRcmdr(c("Population mean != mu0", "Population mean < mu0", 
                                       "Population mean > mu0")), title = gettextRcmdr("Alternative Hypothesis"),
               initialValue = dialog.values$initial.alternative)
  rightFrame <- tkframe(optionsFrame)
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(dialog.values$initial.level)
  confidenceField <- ttkentry(confidenceFrame, width = "6", 
                              textvariable = confidenceLevel)
  muFrame <- tkframe(rightFrame)
  muVariable <- tclVar(dialog.values$initial.mu)
  muField <- ttkentry(muFrame, width = "8", textvariable = muVariable)
  tkgrid(getFrame(xBox), sticky = "nw")
  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
  tkgrid(labelRcmdr(muFrame, text = gettextRcmdr("Null hypothesis: mu = ")), 
         muField, sticky = "w", padx=c(10, 0))
  tkgrid(muFrame, sticky = "w")
  tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: ")), 
         confidenceField, sticky = "w", padx=c(10, 0))
  tkgrid(confidenceFrame, sticky = "w")
  tkgrid(alternativeFrame, rightFrame, sticky = "nw")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
  tkgrid.configure(confidenceField, sticky = "e")
  dialogSuffix()
}

oneWayAnova <- function () {
  Library("multcomp")
  Library("abind")
  defaults <- list(initial.group = NULL, initial.response = NULL, initial.pairwise = 0)
  dialog.values <- getDialog("oneWayAnova", defaults)
  initializeDialog(title = gettextRcmdr("One-Way Analysis of Variance"))
  UpdateModelNumber()
  modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), 
                            sep = ""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
  dataFrame <- tkframe(top)
  groupBox <- variableListBox(dataFrame, Factors(), title = gettextRcmdr("Groups (pick one)"), 
                              initialSelection = varPosn(dialog.values$initial.group, "factor"))
  responseBox <- variableListBox(dataFrame, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"))
  optionsFrame <- tkframe(top)
  pairwiseVariable <- tclVar(dialog.values$initial.pairwise)
  pairwiseCheckBox <- ttkcheckbutton(optionsFrame, variable = pairwiseVariable)
  onOK <- function() {
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)) {
      UpdateModelNumber(-1)
      errorCondition(recall = oneWayAnova, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
                                                             modelValue))
      return()
    }
    if (is.element(modelValue, listAOVModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
        UpdateModelNumber(-1)
        tkdestroy(top)
        oneWayAnova()
        return()
      }
    }
    group <- getSelection(groupBox)
    response <- getSelection(responseBox)
    closeDialog()
    if (length(group) == 0) {
      errorCondition(recall = oneWayAnova, message = gettextRcmdr("You must select a groups factor."))
      return()
    }
    if (length(response) == 0) {
      errorCondition(recall = oneWayAnova, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    command <- paste(modelValue, " <- aov(", response, " ~ ", 
                     group, ", data=", .activeDataSet, ")", sep = "")
    justDoIt(command)
    logger(command)
    doItAndPrint(paste("summary(", modelValue, ")", sep = ""))
    doItAndPrint(paste("with(", .activeDataSet, ", numSummary(",
                       response, ", groups=", group, 
                       ", statistics=c(\"mean\", \"sd\")))", sep = ""))
    activeModel(modelValue)
    putRcmdr("modelWithSubset", FALSE)
    pairwise <- tclvalue(pairwiseVariable)
    putDialog ("oneWayAnova", list (initial.group = group, initial.response = response, initial.pairwise = pairwise))
    if (pairwise == 1) {
      if (eval(parse(text = paste("length(levels(", .activeDataSet, 
                                  "$", group, ")) < 3")))) 
        Message(message = gettextRcmdr("Factor has fewer than 3 levels; pairwise comparisons omitted."), 
                type = "warning")
      else {
        commands <- character(7)
        commands[1] <- paste("local({\n  .Pairs <- glht(", modelValue, 
                             ", linfct = mcp(", group, " = \"Tukey\"))", 
                             sep = "")
        commands[2] <- "  print(summary(.Pairs)) # pairwise tests"
        commands[3] <- "  print(confint(.Pairs)) # confidence intervals"
        commands[4] <- "  print(cld(.Pairs)) # compact letter display"
        commands[5] <- "  old.oma <- par(oma=c(0,5,0,0))"
        commands[6] <- "  plot(confint(.Pairs))"
        commands[7] <- "  par(old.oma)\n})"
        doItAndPrint(paste(commands, collapse="\n"))
      }
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "anova", model = TRUE, reset = "oneWayAnova", apply = "oneWayAnova")
  tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model: ")), 
         model, sticky = "w")
  tkgrid(modelFrame, sticky = "w", columnspan = 2)
  tkgrid(getFrame(groupBox), labelRcmdr(dataFrame, text="  "), getFrame(responseBox), sticky = "nw")
  tkgrid(dataFrame, sticky="w")
  tkgrid(pairwiseCheckBox, labelRcmdr(optionsFrame, text = gettextRcmdr("Pairwise comparisons of means")), 
         sticky = "w")
  tkgrid(optionsFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

multiWayAnova <- function () {
  defaults <- list(initial.group = NULL, initial.response = NULL)
  dialog.values <- getDialog("multiWayAnova", defaults)
  initializeDialog(title = gettextRcmdr("Multi-Way Analysis of Variance"))
  UpdateModelNumber()
  modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), 
                            sep = ""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
  dataFrame <- tkframe(top)
  groupBox <- variableListBox(dataFrame, Factors(), selectmode = "multiple", 
                              title = gettextRcmdr("Factors (pick one or more)"), 
                              initialSelection = varPosn(dialog.values$initial.group, "factor"))
  responseBox <- variableListBox(dataFrame, Numeric(), title = gettextRcmdr("Response Variable (pick one)"), 
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"))
  onOK <- function() {
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)) {
      UpdateModelNumber(-1)
      errorCondition(recall = multiWayAnova, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
                                                               modelValue))
      return()
    }
    if (is.element(modelValue, listAOVModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
        UpdateModelNumber(-1)
        tkdestroy(top)
        multiWayAnova()
        return()
      }
    }
    groups <- getSelection(groupBox)
    response <- getSelection(responseBox)
    putDialog ("multiWayAnova", list (initial.group = groups, initial.response = response))
    closeDialog()
    if (length(groups) == 0) {
      errorCondition(recall = multiWayAnova, message = gettextRcmdr("You must select at least one factor."))
      return()
    }
    if (length(response) == 0) {
      errorCondition(recall = multiWayAnova, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    groups.list <- paste(paste(groups, sep = ""), collapse = ", ")
    doItAndPrint(paste(modelValue, " <- lm(", response, 
                       " ~ ", paste(groups, collapse = "*"), ", data=", 
                       .activeDataSet, ", contrasts=list(", paste(paste(groups, '="contr.Sum"'), collapse=", "), "))", sep = ""))
    doItAndPrint(paste("Anova(", modelValue, ")", sep = ""))
    doItAndPrint(paste("with(", .activeDataSet, ", (tapply(", response, 
                       ", list(", groups.list, "), mean, na.rm=TRUE))) # means", 
                       sep = ""))
    doItAndPrint(paste("with(", .activeDataSet, ", (tapply(", response, 
                       ", list(", groups.list, "), sd, na.rm=TRUE))) # std. deviations", 
                       sep = ""))
    # doItAndPrint(paste("with(", .activeDataSet, ", (tapply(", response, 
    #                    ", list(", groups.list, "), function(x) sum(!is.na(x))))) # counts", 
    #                    sep = ""))
    doItAndPrint(paste("xtabs(~ ", paste(groups, collapse=" + "), ", data=", .activeDataSet, ") # counts", sep=""))
    activeModel(modelValue)
    putRcmdr("modelWithSubset", FALSE)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "Anova", model = TRUE, reset = "multiWayAnova", apply = "multiWayAnova")
  tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model: ")), 
         model, sticky = "w")
  tkgrid(modelFrame, sticky = "w")
  tkgrid(getFrame(groupBox), labelRcmdr(dataFrame, text="  "), getFrame(responseBox), sticky = "nw")
  tkgrid(dataFrame, sticky="w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}
