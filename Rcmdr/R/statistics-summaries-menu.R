# Statistics Menu dialogs

# last modified 2014-08-17 by J. Fox

# Summaries menu

summarizeDataSet <- function(){
	nvar <- length(Variables())
	.activeDataSet <- ActiveDataSet()
	if (nvar > 10){
		response <- RcmdrTkmessageBox(message=sprintf(gettextRcmdr("There are %d variables in the data set %s.\nDo you want to proceed?"), nvar, .activeDataSet),
				icon="question", type="okcancel", default="cancel")
		if ("cancel" == tclvalue(response)) {
			tkfocus(CommanderWindow())
			return()
		}
	}
	doItAndPrint(paste("summary(", .activeDataSet, ")", sep=""))
}

numericalSummaries <- function(){
    Library("abind")
    Library("e1071")
    defaults <- list(initial.x=NULL, initial.mean="1", initial.sd="1", initial.se.mean="0", initial.IQR="1", initial.cv="0",
        initial.quantiles.variable="1", 
        initial.quantiles="0, .25, .5, .75, 1", 
        initial.skewness="0", initial.kurtosis="0", initial.type="2",
        initial.group=NULL, initial.tab=0)
    dialog.values <- getDialog("numericalSummaries", defaults)
    initial.group <- dialog.values$initial.group
    initializeDialog(title=gettextRcmdr("Numerical Summaries"), use.tabs=TRUE, tabs=c("dataTab", "statisticsTab"))
    xBox <- variableListBox(dataTab, Numeric(), selectmode="multiple", title=gettextRcmdr("Variables (pick one or more)"),
        initialSelection=varPosn(dialog.values$initial.x, "numeric"))
    checkBoxes(window = statisticsTab, frame="checkBoxFrame", boxes=c("mean", "sd", "se.mean", "IQR", "cv"), 
        initialValues=c(dialog.values$initial.mean, dialog.values$initial.sd, dialog.values$initial.se.mean, dialog.values$initial.IQR, dialog.values$initial.cv), 
        labels=gettextRcmdr(c("Mean", "Standard Deviation", "Standard Error of Mean", "Interquartile Range", "Coefficient of Variation")))
    skFrame <- tkframe(statisticsTab)
    checkBoxes(window = skFrame, frame="skCheckBoxFrame", boxes=c("skewness", "kurtosis"), 
        initialValues=c(dialog.values$initial.skewness, dialog.values$initial.kurtosis), 
        labels=gettextRcmdr(c("Skewness", "Kurtosis")))
    radioButtons(window = skFrame, name="typeButtons", buttons=c("b1", "b2", "b3"), values=c("1", "2", "3"), 
        initialValue=dialog.values$initial.type,
        labels=gettextRcmdr(c("Type 1", "Type 2", "Type 3")))
    quantilesVariable <- tclVar(dialog.values$initial.quantiles.variable)
    quantilesFrame <- tkframe(statisticsTab)
    quantilesCheckBox <- tkcheckbutton(quantilesFrame, variable=quantilesVariable, 
        text=gettextRcmdr("Quantiles:"))
    quantiles <- tclVar(dialog.values$initial.quantiles)
    quantilesEntry <- ttkentry(quantilesFrame, width="20", textvariable=quantiles)
    groupsBox(recall=numericalSummaries, label=gettextRcmdr("Summarize by:"), 
        initialLabel=if (is.null(initial.group)) gettextRcmdr("Summarize by groups") 
        else paste(gettextRcmdr("Summarize by:"), initial.group), 
        initialGroup=initial.group, window = dataTab)
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        quants <- tclvalue(quantiles)
        meanVar <- tclvalue(meanVariable)
        sdVar <- tclvalue(sdVariable)
        se.meanVar <- tclvalue(se.meanVariable)
        IQRVar <- tclvalue(IQRVariable)
        cvVar <- tclvalue(cvVariable)
        quantsVar <- tclvalue(quantilesVariable)
        skewnessVar <- tclvalue(skewnessVariable)
        kurtosisVar <- tclvalue(kurtosisVariable)
        typeVar <- tclvalue(typeButtonsVariable)
        putDialog("numericalSummaries", list(
            initial.x=x, initial.mean=meanVar, initial.sd=sdVar, initial.se.mean=se.meanVar, initial.IQR=IQRVar, initial.cv=cvVar,
            initial.quantiles.variable=quantsVar, initial.quantiles=quants,
            initial.skewness=skewnessVar, initial.kurtosis=kurtosisVar, initial.type=typeVar,
            initial.group=if (.groups != FALSE) .groups else NULL, initial.tab=tab
        ))      
        if (length(x) == 0){
            errorCondition(recall=numericalSummaries, message=gettextRcmdr("You must select a variable."))
            return()
        }
        closeDialog()
        quants <- paste("c(", gsub(",+", ",", gsub(" ", ",", quants)), ")", sep="")
        .activeDataSet <- ActiveDataSet()
        vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
        else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
        vars <- paste(.activeDataSet, "[,", vars, "]", sep="")
        stats <- paste("c(",
            paste(c('"mean"', '"sd"', '"se(mean)"', '"IQR"', '"quantiles"', '"cv"', '"skewness"', '"kurtosis"')
                [c(meanVar, sdVar, se.meanVar, IQRVar, quantsVar, cvVar, skewnessVar, kurtosisVar) == 1], 
                collapse=", "), ")", sep="")
        if (stats == "c()"){
            errorCondition(recall=numericalSummaries, message=gettextRcmdr("No statistics selected."))
            return()
        }
        type.text <- if (skewnessVar == 1 || kurtosisVar == 1) paste(', type="', typeVar, '"', sep="") else ""
        command <- if (.groups != FALSE) {
            grps <- paste(.activeDataSet, "$", .groups, sep="")
            paste("numSummary(", vars, ", groups=", grps, ", statistics=", stats, 
                ", quantiles=", quants, type.text, ")", sep="")
        }
        else  paste("numSummary(", vars, ", statistics=", stats, 
            ", quantiles=", quants, type.text, ")", sep="")
        doItAndPrint(command) 
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="numSummary", reset="numericalSummaries", apply ="numericalSummaries")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(checkBoxFrame, sticky="w")
    tkgrid(skCheckBoxFrame, typeButtonsFrame, sticky="nw")
    tkgrid(skFrame, sticky="w")
    tkgrid(quantilesCheckBox, quantilesEntry, sticky="w")
    tkgrid(quantilesFrame)
    tkgrid(groupsFrame, sticky = "w", padx=6)
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("dataTab", "statisticsTab"), 
        tab.names=c("Data", "Statistics"))
}

frequencyDistribution <- function () {
  defaults <- list (initial.x = NULL, initial.goodnessOfFit = "0")
  dialog.values <- getDialog ("frequencyDistribution", defaults)
  initializeDialog(title = gettextRcmdr("Frequency Distributions"))
  xBox <- variableListBox(top, Factors(), selectmode = "multiple", 
                          title = gettextRcmdr("Variables (pick one or more)"),
                          initialSelection = varPosn (dialog.values$initial.x, "factor"))
  optionsFrame <- tkframe(top)
  goodnessOfFitVariable <- tclVar(dialog.values$initial.goodnessOfFit)
  goodnessOfFitCheckBox <- ttkcheckbutton(optionsFrame, variable = goodnessOfFitVariable)
  onOK <- function() {
    x <- getSelection(xBox)
    if (length(x) == 0) {
      errorCondition(recall = frequencyDistribution, message = gettextRcmdr("You must select a variable."))
      return()
    }
    goodnessOfFit <- tclvalue(goodnessOfFitVariable)
    putDialog ("frequencyDistribution", list (initial.x = x, initial.goodnessOfFit = goodnessOfFit))
    if (length(x) > 1 && goodnessOfFit == "1") {
      errorCondition(recall = frequencyDistribution, message = gettextRcmdr("Goodness-of-fit test not available when more than one variable is selected."))
      return()
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    for (variable in x) {
      command <- paste("table(", variable, ")", sep = "")
      command <- paste("local({\n  .Table <- with(", .activeDataSet, ", ", command, ")", sep="")
      command <- paste(command, '\n  cat("\\ncounts:\\n")', sep="")
      command <- paste(command, "\n  print(.Table)", sep="")
      command <- paste(command, '\n  cat("\\npercentages:\\n")', sep="")
      command <- paste(command, "\n  print(round(100*.Table/sum(.Table), 2))", sep="")
      if (goodnessOfFit != 1) {
        command <- paste(command, "\n})", sep="")
        doItAndPrint(command)
      }
    }
    env <- environment()
    if (goodnessOfFit == 1) {
      initializeDialog(subwin, title = gettextRcmdr("Goodness-of-Fit Test"))
      hypothesisFrame <- tkframe(subwin)
      levs <- eval(parse(text = paste("levels(", .activeDataSet, 
                                      "$", x, ")", sep = "")))
      n.levs <- length(levs)
      assign(".entry.1", tclVar(paste("1/", n.levs, sep = "")), 
             envir = env)
      make.entries <- "labelRcmdr(hypothesisFrame, text='Hypothesized probabilities:   ')"
      make.lev.names <- "labelRcmdr(hypothesisFrame, text='Factor levels:')"
      for (i in 1:n.levs) {
        entry.varname <- paste(".entry.", i, sep = "")
        assign(entry.varname, tclVar(paste("1/", n.levs, 
                                           sep = "")), envir = env)
        make.entries <- paste(make.entries, ", ", "ttkentry(hypothesisFrame, width='5', textvariable=", 
                              entry.varname, ")", sep = "")
        make.lev.names <- paste(make.lev.names, ", labelRcmdr(hypothesisFrame, text='", 
                                levs[i], "')", sep = "")
      }
      eval(parse(text = paste("tkgrid(", make.lev.names, 
                              ", sticky='w')", sep = "")), envir = env)
      eval(parse(text = paste("tkgrid(", make.entries, 
                              ", stick='w')", sep = "")), envir = env)
      tkgrid(hypothesisFrame, sticky = "w")
      onOKsub <- function() {
        probs <- rep(NA, n.levs)
        for (i in 1:n.levs) {
          entry.varname <- paste(".entry.", i, sep = "")
          res <- try(entry <- eval(parse(text = eval(parse(text = paste("tclvalue(", 
                                                                        entry.varname, ")", sep = "")), envir = env))), 
                     silent = TRUE)
          if (class(res) == "try-error") {
            errorCondition(subwin, message = gettextRcmdr("Invalid entry."))
            return()
          }
          if (length(entry) == 0) {
            errorCondition(subwin, message = gettextRcmdr("Missing entry."))
            return()
          }
          opts <- options(warn = -1)
          probs[i] <- as.numeric(entry)
          options(opts)
        }
        probs <- na.omit(probs)
        if (length(probs) != n.levs) {
          errorCondition(subwin, message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number levels (%d)."), 
                                                   length(probs), n.levs))
          return()
        }
        if (any(probs < 0)) {
          errorCondition(subwin, message = gettextRcmdr("Negative probabilities not allowed."))
          return()
        }
        if (abs(sum(probs) - 1) > 0.001) {
          Message(message = gettextRcmdr("Probabilities rescaled to sum to 1."), 
                  type = "warning")
          probs <- probs/sum(probs)
        }
        closeDialog(subwin)
        command <- paste(command, "\n  .Probs <- c(", paste(probs, collapse = ","), ")", sep = "")
        command <- paste(command, "\n  chisq.test(.Table, p=.Probs)\n})")
        doItAndPrint(command)
      }
      subOKCancelHelp(subwin)
      tkgrid(subButtonsFrame, sticky = "w")
      dialogSuffix(subwin, onOK = onOKsub, focus = subwin, force.wait=TRUE)
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "table", reset = "frequencyDistribution", apply="frequencyDistribution")
  tkgrid(getFrame(xBox), sticky = "nw")
  tkgrid(goodnessOfFitCheckBox, 
         labelRcmdr(optionsFrame, text = gettextRcmdr("Chi-square goodness-of-fit test (for one variable only)")), 
         sticky = "w")
  tkgrid(optionsFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

statisticsTable <- function () {
  defaults <- list (initial.group=NULL, initial.response=NULL, initial.statistic="mean", initial.other = "")
  dialog.values <- getDialog ("statisticsTable", defaults)
  initializeDialog(title = gettextRcmdr("Table of Statistics"))
  variablesFrame <- tkframe(top)
  groupBox <- variableListBox(variablesFrame, Factors(), selectmode = "multiple", 
                              title = gettextRcmdr("Factors (pick one or more)"), 
                              initialSelection = varPosn(dialog.values$initial.group,"factor"))
  responseBox <- variableListBox(variablesFrame, Numeric(), selectmode = "multiple", 
                                 initialSelection = varPosn(dialog.values$initial.response, "numeric"),
                                 title = gettextRcmdr("Response variables (pick one or more)"))
  statFrame <- tkframe(top)
  radioButtons(statFrame, name = "statistic", buttons = c("mean", "median", "sd", "IQR", "other"), 
               labels = gettextRcmdr(c("Mean", "Median", "Standard deviation", "Interquartile range", "Other (specify)")), 
               initialValue = dialog.values$initial.statistic, 
               title = gettextRcmdr("Statistic"))
  otherVariable <- tclVar(dialog.values$initial.other)
  otherEntry <- ttkentry(statFrame, width = "20", textvariable = otherVariable)
  tkgrid(statisticFrame, labelRcmdr(statFrame, text ="  "), otherEntry, sticky = "sw")
  onOK <- function() {
    groups <- getSelection(groupBox)
    if (0 == length(groups)) {
      errorCondition(recall = statisticsTable, message = gettextRcmdr("No factors selected."))
      return()
    }
    responses <- getSelection(responseBox)
    if (0 == length(responses)) {
      errorCondition(recall = statisticsTable, message = gettextRcmdr("You must select a response variable."))
      return()
    }
    stat <- statistic <- tclvalue(statisticVariable)
    if (statistic == "other") 
      statistic <- tclvalue(otherVariable)
    putDialog ("statisticsTable", list(initial.group=groups, initial.response=responses, 
                                       initial.statistic=stat, initial.other = if(stat == "other") statistic else ""))  
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    groups.list <- paste(paste(groups, sep = ""), collapse = ", ")
    for (response in responses) {
      if (length(responses) > 1) 
        doItAndPrint(paste("# Table for ", response, 
                           ":", sep = ""))
      doItAndPrint(paste("with(", .activeDataSet, ", tapply(",  
                         response, ", list(", groups.list, "), ", statistic, 
                         ", na.rm=TRUE))", sep = ""))
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "tapply", reset="statisticsTable", apply="statisticsTable")
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text = "    "), 
         getFrame(responseBox), sticky = "nw")
  tkgrid(variablesFrame, sticky = "w")
  tkgrid(statFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix(focus = otherEntry)
}

correlationMatrix <- function (){
    defaults <- list (initial.x = NULL, initial.correlations = "Pearson", initial.pvaluesVar="0", initial.use="complete")  
    dialog.values <- getDialog ("correlationMatrix", defaults)
    initializeDialog(title = gettextRcmdr("Correlation Matrix"))
    xBox <- variableListBox(top, Numeric(), selectmode = "multiple", 
                            title = gettextRcmdr("Variables (pick two or more)"),
                            initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    radioButtons(name = "correlations", buttons = c("pearson", 
                                                    "spearman", "partial"), values = c("Pearson", "Spearman", 
                                                                                       "partial"), labels = gettextRcmdr(c("Pearson product-moment", 
                                                                                                                           "Spearman rank-order", "Partial")), title = gettextRcmdr("Type of Correlations"),
                 initialValue = dialog.values$initial.correlations)
    radioButtons(name = "use", buttons = c("complete", "pairwise.complete"), 
                 labels = gettextRcmdr(c("Complete observations", "Pairwise-complete observations")), title = gettextRcmdr("Observations to Use"),
                 initialValue = dialog.values$initial.use)
    pvaluesFrame <- tkframe(top)
    pvaluesVar <- tclVar(dialog.values$initial.pvaluesVar)
    pvaluesCheckbox <- ttkcheckbutton(pvaluesFrame, variable = pvaluesVar, text = gettextRcmdr("Pairwise p-values"))
    onOK <- function() {
        correlations <- tclvalue(correlationsVariable)
        use <- tclvalue(useVariable)
        x <- getSelection(xBox)
        pvalues <- tclvalue(pvaluesVar)
        if (2 > length(x)) {
            errorCondition(recall = correlationMatrix, message = gettextRcmdr("Fewer than 2 variables selected."))
            return()
        }
        if ((correlations == "partial") && (3 > length(x))) {
            errorCondition(recall = correlationMatrix, message = gettextRcmdr("Fewer than 3 variables selected\nfor partial correlations."))
            return()
        }
        closeDialog()
        putDialog ("correlationMatrix", list (initial.x=x, initial.correlations=correlations, 
                                              initial.pvaluesVar=pvalues, initial.use=use))
        x <- paste("\"", x, "\"", sep = "")
        .activeDataSet <- ActiveDataSet()
        if (correlations == "Pearson") {
            if (pvalues == 0) {
                doItAndPrint(paste("cor(", .activeDataSet, "[,c(", 
                                   paste(x, collapse = ","), ")], use=\"", use, "\")", 
                                   sep = ""))
            }
            else {
                Library("Hmisc")
                doItAndPrint(paste("rcorr.adjust(", .activeDataSet, 
                                   "[,c(", paste(x, collapse = ","), ")], type=\"pearson\", use=\"", use, "\")", 
                                   sep = ""))
            }
        }
        else if (correlations == "Spearman") {
            if (pvalues == 0) {
                doItAndPrint(paste("cor(", .activeDataSet, "[,c(", 
                                   paste(x, collapse = ","), ")], method=\"spearman\", use=\"", use, "\")", 
                                   sep = ""))
            }
            else {
                Library("Hmisc")
                doItAndPrint(paste("rcorr.adjust(", .activeDataSet, 
                                   "[,c(", paste(x, collapse = ","), ")], type=\"spearman\", use=\"", use, "\")", 
                                   sep = ""))
            }
        }
        else if (pvalues == 0){
            doItAndPrint(paste("partial.cor(", .activeDataSet, 
                               "[,c(", paste(x, collapse = ","), ")], use=\"", use, "\")", 
                               sep = ""))
        }
        else {
            doItAndPrint(paste("partial.cor(", .activeDataSet, 
                               "[,c(", paste(x, collapse = ","), ")], tests=TRUE, use=\"", use, "\")", 
                               sep = ""))
        }
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "rcorr.adjust", reset="correlationMatrix", apply="correlationMatrix")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(correlationsFrame, sticky = "w")
    tkgrid(useFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=""))
    tkgrid(pvaluesCheckbox, sticky="w") 
    tkgrid(pvaluesFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix()
}

# the following dialog contributed by Stefano Calza, modified by J. Fox

correlationTest <- function(){
  defaults <- list(initial.x=NULL,initial.correlations="pearson",initial.alternative ="two.sided")
  dialog.values <- getDialog("correlationTest", defaults)
  initializeDialog(title=gettextRcmdr("Correlation Test"))
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettextRcmdr("Variables (pick two)"),initialSelection=varPosn(dialog.values$initial.x, "numeric"))
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="correlations", buttons=c("pearson", "spearman", "kendall"),
               labels=gettextRcmdr(c("Pearson product-moment", "Spearman rank-order", "Kendall's tau")),
               initialValue=dialog.values$initial.correlations, 
               title=gettextRcmdr("Type of Correlation"))
  radioButtons(optionsFrame, name="alternative", buttons=c("two.sided", "less", "greater"), 
               values=c("two.sided", "less", "greater"),
               initialValue=dialog.values$initial.alternative, 
               labels=gettextRcmdr(c("Two-sided", "Correlation < 0", "Correlation > 0")), 
               title=gettextRcmdr("Alternative Hypothesis"))  
  onOK <- function(){
    alternative <- as.character(tclvalue(alternativeVariable))
    correlations <- as.character(tclvalue(correlationsVariable))
    x <- getSelection(xBox)
    putDialog("correlationTest", list(initial.alternative=alternative, initial.correlations=correlations, initial.x=x))
    if (2 > length(x)) {
      errorCondition(recall=correlationTest,
                     message=gettextRcmdr("Fewer than 2 variables selected."))
      return()
    }
    if(2 < length(x)) {
      errorCondition(recall=correlationTest,
                     message=gettextRcmdr("More than 2 variables selected."))
      return()
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    command <- paste("with(", .activeDataSet, ", cor.test(", x[1], ", ", x[2],
                     ', alternative="', alternative, '", method="', correlations, '"))', sep="")
    doItAndPrint(command)  
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cor.test", reset="correlationTest", apply="correlationTest")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(correlationsFrame, labelRcmdr(optionsFrame, text="  "), alternativeFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix()
}

countMissing <- function(){
	command <- paste("sapply(", activeDataSet(), 
			", function(x)(sum(is.na(x)))) # NA counts", sep="")
	doItAndPrint(command)
	invisible(NULL)
}

ShapiroTest <- function () {
  defaults <- list (initial.var = NULL)
  dialog.values <- getDialog ("ShapiroTest", defaults)
  initializeDialog(title = gettextRcmdr("Shapiro-Wilk Test for Normality"))
  variableBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"),
                                 initialSelection = varPosn (dialog.values$initial.var, "numeric"))
  onOK <- function() {
    var <- getSelection(variableBox)
    putDialog ("ShapiroTest", list (initial.var = var))
    if (length(var) == 0) {
      errorCondition(recall = ShapiroTest, message = gettextRcmdr("You must select a variable."))
      return()
    }
    closeDialog()
    doItAndPrint(paste("with(", ActiveDataSet(), ", shapiro.test(", 
                       var, "))", sep = ""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "shapiro.test", reset = "ShapiroTest", apply = "ShapiroTest")
  tkgrid(getFrame(variableBox), sticky = "nw")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}
