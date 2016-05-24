
#=========================================================================================================================================

fncvROCR <- function(){
  varPosListn <- function(vars, var){
    if (is.null(var)) return(NULL)
    if (any(!var %in% vars)) NULL
    else apply(outer(var, vars, "=="), 1, which) - 1
  }

  #require(ROCR)
  #Daniel
  performancelist <- c("acc", "err", "fpr", "fall", "tpr", "rec", "sens", "fnr", "miss", "tnr", "spec", "ppv", 
                       "prec", "npv", "pcfall", "pcmiss", "rpp", "rnp", "phi", "mat", "mi", "chisq", "odds", 
                       "lift", "f", "rch", "auc", "prbe", "cal", "mxe", "rmse", "sar", "ecost", "cost")
  performancelistlong <- c("Accuracy", "Error rate", "False positive rate", "Fallout (fpr)", "True positive rate", "Recall (tpr)", "Sensitivity", "False negative rate", "Miss (fnr)", "True negative rate", "Specificity", "Positive predictive value", 
                       "Precision (ppv)", "Negative predictive value", "Prediction-conditioned fallout", "Prediction-conditioned miss", "Rate of positive predictions", "Rate of negative predictions", "Phi correlation coefficient", "Mattheus correlation coefficient (phi)", "Mutual information", "Chi square test statistic", "Odds ratio", 
                       "Lift value", "Precision-recall F measure", "ROC convex hull", "AUC", "Precision-recall break-even point", "Callibration error", "Mean cross-entropy", "Root-mean-squared error", "Sar", "Expected cost", "Cost")
  defaults <- list(initial.prediction = NULL, initial.label = NULL, initial.ymeasure = performancelistlong[5], initial.xmeasure = performancelistlong[3],
                   initial.colorize = 0, initial.add = 0,
                   initial.printcutoffs = 0, initial.cutoffs = "seq(0,1,by=0.1)",
                   initial.printroc = 0, 
                   initial.costfp = 1, initial.costfn = 1, 
                   initial.calwindowsize = 100, 
                   initial.partialfprstop = 1,
                   initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"), 
                   initial.main=gettextRcmdr("<auto>"),
                   initial.tab=0) # tab
  dialog.values <- getDialog("ROCR", defaults)
  
  initializeDialog(title=gettext("Plot ROC curve", domain="R-RcmdrPlugin.ROC"), use.tabs=TRUE) # tab

  #Daniel
  .factors <- Factors()
  .numeric <- Numeric()
  predictionBox <- variableListBox(dataTab, .numeric, title=gettext("Predictions variable (pick one)", domain="R-RcmdrPlugin.ROC"),# tab
                                   initialSelection=varPosn(dialog.values$initial.prediction, "numeric"))
  labelBox <- variableListBox(dataTab, .numeric, title=gettext("Labels variable (pick one)", domain="R-RcmdrPlugin.ROC"),
                              initialSelection=varPosn(dialog.values$initial.label, "numeric"))

  optionsParFrame <- tkframe(optionsTab)# tab
  parFrame <- ttklabelframe(optionsParFrame, text=gettext("Plot Labels and Points", domain="R-RcmdrPlugin.ROC"))# tab
  performanceFrame <- ttklabelframe(optionsParFrame, text=gettext("Performance measures", domain="R-RcmdrPlugin.ROC"))# tab
  #performanceoptFrame <- ttklabelframe(optionsParFrame, text=gettext("Performance options", domain="R-RcmdrPlugin.ROC"))# tab
  
  costfpVar <- tclVar(dialog.values$initial.costfp) # tab
  costfpEntry <- ttkentry(performanceFrame, width = "25", textvariable = costfpVar)# tab
  costfnVar <- tclVar(dialog.values$initial.costfn) # tab
  costfnEntry <- ttkentry(performanceFrame, width = "25", textvariable = costfnVar)# tab
  calwindowsizeVar <- tclVar(dialog.values$initial.calwindowsize) # tab
  calwindowsizeEntry <- ttkentry(performanceFrame, width = "25", textvariable = calwindowsizeVar)# tab
  fprstopVar <- tclVar(dialog.values$initial.partialfprstop) # tab
  fprstopEntry <- ttkentry(performanceFrame, width = "25", textvariable = fprstopVar)# tab
  
  checkBoxes(window = optionsParFrame, frame = "optionsFrame",# tab
             boxes = c("printroc", "colorize", "add", "printcutoffs"), initialValues = c(
               dialog.values$initial.printroc, dialog.values$initial.colorize, dialog.values$initial.add, dialog.values$initial.printcutoffs),labels = gettextRcmdr(c(
                 "Print performance object", "Colorize according to cutoff", "Add curve to existing plot","Print cutoffs")), title = gettext("Plot Options", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  
  ymeasureBox <- variableListBox(performanceFrame, performancelistlong, title=gettext("Performance measure (y) (pick one)", domain="R-RcmdrPlugin.ROC"),# tab
                                 initialSelection=varPosListn(performancelistlong, dialog.values$initial.ymeasure))
  xmeasureBox <- variableListBox(performanceFrame, performancelistlong, title=gettext("Performance measure (x) (pick one or none)", domain="R-RcmdrPlugin.ROC"),
                                 initialSelection=varPosListn(performancelistlong, dialog.values$initial.xmeasure))
  
cutoffsVar <- tclVar(dialog.values$initial.cutoffs) # tab
cutoffsEntry <- ttkentry(optionsFrame, width = "25", textvariable = cutoffsVar)# tab
cutoffsScroll <- ttkscrollbar(optionsFrame, orient = "horizontal",
                           command = function(...) tkxview(cutoffsEntry, ...))
tkconfigure(cutoffsEntry, xscrollcommand = function(...) tkset(cutoffsScroll,
                                                            ...))
tkbind(cutoffsEntry, "<FocusIn>", function() tkselection.clear(cutoffsEntry))
tkgrid(labelRcmdr(optionsFrame, text = gettext("Print cutoffs at", domain="R-RcmdrPlugin.ROC")), cutoffsEntry, sticky = "ew", padx=6)
tkgrid(labelRcmdr(optionsFrame, text =""), cutoffsScroll, sticky = "ew", padx=6)

  
xlabVar <- tclVar(dialog.values$initial.xlab) # tab
ylabVar <- tclVar(dialog.values$initial.ylab)
mainVar <- tclVar(dialog.values$initial.main)
xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                           command = function(...) tkxview(xlabEntry, ...))
tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
                                                            ...))
tkbind(xlabEntry, "<FocusIn>", function() tkselection.clear(xlabEntry))
tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                           command = function(...) tkxview(ylabEntry, ...))
tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
                                                            ...))
tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
tkgrid(labelRcmdr(parFrame, text=""), ylabScroll, sticky = "ew", padx=6)
mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                           command = function(...) tkxview(mainEntry, ...))
tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
                                                            ...))
tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)

  onOK <- function(){
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1 # tab
    #Daniel
    prediction <- getSelection(predictionBox)
    label <- getSelection(labelBox)    
    ymeasure <- getSelection(ymeasureBox)
    xmeasure <- getSelection(xmeasureBox)  
    
    colorize <- as.character("1" == tclvalue(colorizeVariable))
    add <- as.character("1" == tclvalue(addVariable))
    printroc <- as.character("1" == tclvalue(printrocVariable))

    costfp = as.numeric(as.character(tclvalue(costfpVar)))
    costfn = as.numeric(as.character(tclvalue(costfnVar)))
    calwindowsize = as.numeric(as.character(tclvalue(calwindowsizeVar)))
    fprstop = as.numeric(as.character(tclvalue(fprstopVar)))
    
    printcutoffsat <- if ("0" == tclvalue(printcutoffsVariable))
      ""
    else paste(", print.cutoffs.at=", tclvalue(cutoffsVar), sep = "")    

    xlab <- trim.blanks(tclvalue(xlabVar))
    xlab <- if (xlab == gettextRcmdr("<auto>"))
      ""
    else paste(", xlab=\"", xlab, "\"", sep = "")
    ylab <- trim.blanks(tclvalue(ylabVar))
    ylab <- if (ylab == gettextRcmdr("<auto>"))
      ""
    else paste(", ylab=\"", ylab, "\"", sep = "")
    main <- trim.blanks(tclvalue(mainVar))
    main <- if (main == gettextRcmdr("<auto>"))
      ""
    else paste(", main=\"", main, "\"", sep = "")
   

putDialog ("ROCR", list(initial.prediction = prediction, initial.label = label, initial.ymeasure = ymeasure, initial.xmeasure = xmeasure,
                        initial.colorize = tclvalue(colorizeVariable), initial.add = tclvalue(addVariable),
                        initial.printcutoffs = tclvalue(printcutoffsVariable), initial.cutoffs = tclvalue(cutoffsVar),
                        initial.printroc = tclvalue(printrocVariable), 
                        initial.costfp = as.numeric(as.character(tclvalue(costfpVar))),
                        initial.costfn = as.numeric(as.character(tclvalue(costfnVar))),
                        initial.calwindowsize = as.numeric(as.character(tclvalue(calwindowsizeVar))),
                        initial.partialfprstop = as.numeric(as.character(tclvalue(fprstopVar))),
                        initial.xlab=tclvalue(xlabVar), initial.ylab=tclvalue(ylabVar), 
                        initial.main=tclvalue(mainVar),
                        initial.tab=tab)) # tab
closeDialog()
    
   if (0 == length(prediction)) {
      errorCondition(recall=fncvROCR, message=gettext("You must select a prediction variable.", domain="R-RcmdrPlugin.ROC"))
      return()
    }
    if (0 == length(label)) {
      errorCondition(recall=fncvROCR, message=gettext("No labels variables selected.", domain="R-RcmdrPlugin.ROC"))
      return()
    }
   if (0 == length(ymeasure)) {
     errorCondition(recall=fncvROCR, message=gettext("You must select a performance measure (y) variable.", domain="R-RcmdrPlugin.ROC"))
     return()
   }
   
    if (0 != length(xmeasure)) {
      if (ymeasure == xmeasure) {
       errorCondition(recall=fncvROCR, message=gettext("The performance measures, x and y should be different.", domain="R-RcmdrPlugin.ROC"))
       return()
     }
    }

    .activeDataSet <- ActiveDataSet()
    

    #Daniel
    command <- paste("pred <- prediction(", .activeDataSet, "$", prediction, ", ", 
                     .activeDataSet, "$", label, ")", sep = "")
    doItAndPrint(command)

    ymeasure <- performancelist[which(performancelistlong == ymeasure)]

    if (ymeasure == "auc") {
      .partialfprstop <- paste(", fpr.stop=", fprstop, sep = "")
    } else {
      .partialfprstop <- ""
    }
    if (ymeasure == "cal") {
      .calwindowsize <- paste(", window.size=", calwindowsize, sep = "")
    } else {
      .calwindowsize <- ""
    }
    if (ymeasure == "cost") {
      .cost <- paste(", cost.fp=", costfp, ", cost.fn=", costfn, sep = "")
    } else {
      .cost <- ""
    }

    if (0 == length(xmeasure)) {
      command <- paste("perf <- performance(pred, '", ymeasure, "'", .partialfprstop, .calwindowsize, .cost, ")", sep = "")
      doItAndPrint(command)
    } else {
      command <- paste("perf <- performance(pred, '", ymeasure, "', '", 
                       performancelist[which(performancelistlong == xmeasure)], "'", .partialfprstop, .calwindowsize, .cost, ")", sep = "")
      doItAndPrint(command)      
    }

    if (printroc == "TRUE") {
      command <- paste("perf", sep = "")
      doItAndPrint(command)
    }

    command <- paste("plot(perf, colorize=", colorize, ", add=", add, printcutoffsat, xlab, ylab, main, ")", sep = "")
    doItAndPrint(command)

    command <- paste("remove(perf)", sep = "")
    doItAndPrint(command)
    command <- paste("remove(pred)", sep = "")
    doItAndPrint(command)

    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="performance", reset = "fncvROCR", apply="fncvROCR")

tkgrid(getFrame(predictionBox), getFrame(labelBox), sticky = "nw", padx=6, pady=c(6, 0))
tkgrid(getFrame(ymeasureBox), getFrame(xmeasureBox), sticky="nw", padx=6, pady=c(6, 0))
tkgrid(performanceFrame, sticky = "we", padx=6, pady=c(6, 6))
tkgrid(labelRcmdr(performanceFrame, text = gettext("Partial AUC upt to fpr of", domain="R-RcmdrPlugin.ROC")), fprstopEntry, sticky = "ew", padx=6, pady=c(6, 0))
tkgrid(labelRcmdr(performanceFrame, text = gettext("Calibration error window size", domain="R-RcmdrPlugin.ROC")), calwindowsizeEntry, sticky = "ew", padx=6, pady=c(0, 0))
tkgrid(labelRcmdr(performanceFrame, text = gettext("Cost fp adjustment", domain="R-RcmdrPlugin.ROC")), costfpEntry, sticky = "ew", padx=6, pady=c(0, 0))
tkgrid(labelRcmdr(performanceFrame, text = gettext("Cost fn adjustment", domain="R-RcmdrPlugin.ROC")), costfnEntry, sticky = "ew", padx=6, pady=c(0, 6))
#tkgrid(performanceoptFrame, sticky = "we", padx=6, pady=c(6, 6))
#tkgrid(getFrame(performanceFrame), getFrame(performanceoptFrame), sticky="nw", padx=6, pady=c(6, 6))
tkgrid(optionsParFrame, sticky = "we", padx=6, pady=c(6, 0))
tkgrid(optionsFrame, parFrame, sticky = "nswe", padx=6, pady=6)


tkgrid(ttklabel(dataTab, text=""))
tkgrid(ttklabel(dataTab, text=""))
tkgrid(labelRcmdr(top, text = " "), padx=6)
dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}    

#=========================================================================================================================================
