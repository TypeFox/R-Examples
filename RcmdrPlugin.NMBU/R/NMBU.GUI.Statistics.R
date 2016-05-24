# $Id: NMBU.GUI.Statistics.R 35 2014-01-10 21:17:26Z khliland $

##
## GUI functions for the Statistics menu
##



## GUI tips
#
# Usual code structure:
#    1. Intitialise window and prepare graphical elements
#    2. onOK function contianing actions to perform
#       2.1 Collect values from GUI
#       2.2 Test if combination of values is usable
#       2.3 Perform main calculations, print, update models/data etc.
#    3. Set up GUI.
#       - tkgrid() adds elements. Explicit placement and width/heigth by colum, row, columnspan and rowspan
#       - Frames with graphical elements are safer than direct placement of elements due to version conflicts.
#       - dialogSuffix() defines the final size of the grid used for elements.

############################
# Covariance matrix
covarianceMatrix <- function(){
  initializeDialog(title=gettextRcmdr("Covariance/correlation matrix")) # Window heading
  initial.group <- NULL
  # Prepare selection box and radio buttons
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettextRcmdr("Variables (pick two or more)"))
  wBox <- variableListBox(top, c("-none-", Numeric()), selectmode="single", title=gettextRcmdr("Weights (pick zero or one)"))
  radioButtons(name="covariances", buttons=c("pearson", "spearman"), values=c("Pearson", "Spearman"),
               labels=gettextRcmdr(c("Pearson product-moment", "Spearman rank-order")), title=gettextRcmdr("Type (unweighted)"))
  radioButtons(name="covCor", buttons=c("covariance", "correlation"), values=c("covariance", "correlation"),
               labels=gettextRcmdr(c("Covariance", "Correlation")), title=gettextRcmdr("Output"))
  groupsBox(recall=covarianceMatrix, label=gettextRcmdr("Compute by:"), 
            initialLabel=gettextRcmdr("Compute by groups"), 
            initialGroup=initial.group)
  
  onOK <- function(){ # Actions to perform
    covariances <- tclvalue(covariancesVariable)
    covCor <- tclvalue(covCorVariable)
    x <- getSelection(xBox)
    w <- getSelection(wBox)
    if (2 > length(x)) {
      errorCondition(recall=covarianceMatrix, message=gettextRcmdr("Fewer than 2 variables selected."))
      return()
    }
    closeDialog()
    x <- paste('"', x, '"', sep="")
    .activeDataSet <- ActiveDataSet()
    if(.groups==FALSE){
      if(length(w)>0 && w !="-none-"){
        if(covCor=="correlation"){
          doItAndPrint(paste("cov.wt(", .activeDataSet, "[,c(", paste(x, collapse=","),
                             ")], wt=", .activeDataSet, "$", w ,", cor=TRUE)$cor", sep=""))
        } else {
          doItAndPrint(paste("cov.wt(", .activeDataSet, "[,c(", paste(x, collapse=","),
                             ")], wt=", .activeDataSet, "$", w ,")$cov", sep=""))
        }
      } else {
        type <- ifelse(covCor=="correlation","cor","cov")
        if (covariances == "Pearson"){
          doItAndPrint(paste(type, "(", .activeDataSet, "[,c(", paste(x, collapse=","),
                             ')], use="complete.obs")', sep=""))
        }
        else if (covariances == "Spearman"){
          logger("# Spearman rank-order covariances")
          doItAndPrint(paste(type, "(", .activeDataSet, "[,c(", paste(x, collapse=","),
                             ')], use="complete.obs", method="spearman")', sep=""))
        }
      }
    } else {
      eval(parse(text=paste(".levels <- levels(", .activeDataSet, "$", .groups, ")", sep="")))
      for(i in 1:length(.levels)){
        if(length(w)>0 && w !="-none-"){
          if(covCor=="correlation"){
            doItAndPrint(paste("cov.wt(", .activeDataSet, "[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"',c(", paste(x, collapse=","),
                               ")], wt=", .activeDataSet, "$", w ,"[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"'], cor=TRUE)$cor", sep=""))
          } else {
            doItAndPrint(paste("cov.wt(", .activeDataSet, "[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"',c(", paste(x, collapse=","),
                               ")], wt=", .activeDataSet, "$", w ,"[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"'])$cov", sep=""))
          }
        } else {
          type <- ifelse(covCor=="correlation","cor","cov")
          if (covariances == "Pearson"){
            doItAndPrint(paste(type, "(", .activeDataSet, "[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"',c(", paste(x, collapse=","),
                               ')], use="complete.obs")', sep=""))
          }
          else if (covariances == "Spearman"){
            logger("# Spearman rank-order covariances")
            doItAndPrint(paste(type, "(", .activeDataSet, "[", .activeDataSet, "$", .groups, "=='", .levels[i] ,"',c(", paste(x, collapse=","),
                               ')], use="complete.obs", method="spearman")', sep=""))
          }
        }
      }
    }
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="cov")
  tkgrid(getFrame(xBox), sticky="nw", row=1, column=1)
  tkgrid(getFrame(wBox), sticky="nw", row=1, column=2)
  tkgrid(groupsFrame, sticky = "w", row=2, column=1)
  tkgrid(covariancesFrame, sticky="w", row=3, column=1)
  tkgrid(covCorFrame, sticky="w", row=3, column=2)
  tkgrid(buttonsFrame, sticky="w", row=4, column=1, columnspan=2)
  dialogSuffix(rows=4, columns=2)
}

####################################
# Proportion testing (two proportions) (without data)
twoProportionTest <- function(){
  initializeDialog(title=gettextRcmdr("Two-Sample Proportion Test"))
  onOK <- function(){ # Actions to perform
    succ1 <- tclvalue(successLevel1)
    succ2 <- tclvalue(successLevel2)
    if(trim.blanks(succ1) == gettextRcmdr("") || trim.blanks(succ2) == gettextRcmdr("")){
      errorCondition(recall=twoProportionTest, message=gettextRcmdr("Please specify the number of successes for both groups."))
      return()
    }
    fail1 <- tclvalue(failureLevel1)
    fail2 <- tclvalue(failureLevel2)
    if(trim.blanks(fail1) == gettextRcmdr("") || trim.blanks(fail2) == gettextRcmdr("")){
      errorCondition(recall=twoProportionTest, message=gettextRcmdr("Please specify the number of failures for both groups."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    test    <- as.character(tclvalue(testVariable))
    testVer <- as.character(tclvalue(testVerVariable))
    closeDialog()
    command <- paste(".Table <- cbind(successes=c(", succ1, ",", succ2, "), failures=c(", fail1, ",", fail2, "))", sep="")
    #    logger(command)
    #    assign(".Table", justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(command)
    doItAndPrint(".Table")
    if(testVer=="pooled"){
      pooled <- "TRUE)"
    } else {
      pooled <- "FALSE)"
    }
    if (test == "normal") doItAndPrint(paste("prop.test.ordinary(.Table, alternative='", 
                                             alternative, "', conf.level=", level, ", correct=FALSE, pooled=", pooled, sep=""))
    else doItAndPrint(paste("prop.test.ordinary(.Table, alternative='", 
                            alternative, "', conf.level=", level, ", correct=TRUE, pooled=", pooled, sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  successFrame  <- tkframe(top)
  successLevel1 <- tclVar("?")
  successLevel2 <- tclVar("?")
  successField1 <- ttkentry(successFrame, width="6", textvariable=successLevel1)
  successField2 <- ttkentry(successFrame, width="6", textvariable=successLevel2)
  failureFrame  <- tkframe(top)
  failureLevel1 <- tclVar("?")
  failureLevel2 <- tclVar("?")
  failureField1 <- ttkentry(failureFrame, width="6", textvariable=failureLevel1)
  failureField2 <- ttkentry(failureFrame, width="6", textvariable=failureLevel2)
  tkgrid(labelRcmdr(successFrame, text=gettextRcmdr("# of successes:"), fg="blue"), successField1, successField2, sticky="nw")
  tkgrid(labelRcmdr(failureFrame, text=gettextRcmdr("# of failures:"), fg="blue"), failureField1, failureField2, sticky="nw")
  tkgrid(successFrame, sticky="nw")
  tkgrid(failureFrame, sticky="nw")
  OKCancelHelp(helpSubject="prop.test")
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < 0", "Difference > 0")), title=gettextRcmdr("Alternative Hypothesis"))
  rightFrame <- tkframe(top)
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  radioButtons(name="test", buttons=c("normal", "corrected"),
               labels=gettextRcmdr(c("Normal approximation", "Normal approximation with\ncontinuity correction")), 
               title=gettextRcmdr("Type of Test"))
  radioButtons(name="testVer", buttons=c("pooled", "individual"), 
               labels=gettextRcmdr(c("Pooled sample proportion", "Individual sample proportions")), 
               title=gettextRcmdr("Approximation version"))
  tkgrid(labelRcmdr(rightFrame, text=""))
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level: "), fg="blue"), confidenceField, sticky="nw")
  tkgrid(confidenceFrame, sticky="nw")
  tkgrid(alternativeFrame, rightFrame, sticky="nw")
  tkgrid(testFrame, testVerFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="nw")
  tkgrid.configure(confidenceField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}

####################################
# Proportion testing (without data)
proportionTest <- function(){
  initializeDialog(title=gettextRcmdr("Single-Sample Proportion Test"))
  onOK <- function(){ # Actions to perform
    succ <- tclvalue(successLevel)
    if(trim.blanks(succ) == gettextRcmdr("")){
      errorCondition(recall=proportionTest, message=gettextRcmdr("Please specify the number of successes."))
      return()
    }
    fail <- tclvalue(failureLevel)
    if(trim.blanks(fail) == gettextRcmdr("")){
      errorCondition(recall=proportionTest, message=gettextRcmdr("Please specify the number of failures."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    test    <- as.character(tclvalue(testVariable))
    testVer <- as.character(tclvalue(testVerVariable))
    p <- tclvalue(pVariable)
    closeDialog()
    command <- paste(".Table <- cbind(successes=", succ, ", failures=", fail, ")", sep="")
    #    logger(command)
    #    assign(".Table", justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(command)
    doItAndPrint(".Table")
    if(testVer=="ordinaryVer"){
      the.prop <- "prop.test.ordinary"
    } else {
      the.prop <- "prop.test"
    }
    if (test == "normal") doItAndPrint(paste(the.prop, "(.Table, alternative='", 
                                             alternative, "', p=", p, ", conf.level=", level, ", correct=FALSE)", sep=""))
    else if (test == "corrected") doItAndPrint(paste(the.prop, "(.Table, alternative='", 
                                                     alternative, "', p=", p, ", conf.level=", level, ", correct=TRUE)", sep=""))
    else doItAndPrint(paste("binom.test(.Table, alternative='", 
                            alternative, "', p=", p, ", conf.level=", level, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  successFrame <- tkframe(top)
  successLevel <- tclVar("?")
  successField <- ttkentry(successFrame, width="6", textvariable=successLevel)
  failureFrame <- tkframe(top)
  failureLevel <- tclVar("?")
  failureField <- ttkentry(failureFrame, width="6", textvariable=failureLevel)
  tkgrid(labelRcmdr(successFrame, text=gettextRcmdr("# of successes:"), fg="blue"), successField, sticky="nw")
  tkgrid(labelRcmdr(failureFrame, text=gettextRcmdr("# of failures:"), fg="blue"), failureField, sticky="nw")
  tkgrid(successFrame, sticky="nw")
  tkgrid(failureFrame, sticky="nw")
  OKCancelHelp(helpSubject="prop.test")
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Population proportion != p0", "Population proportion < p0", "Population proportion > p0")), title=gettextRcmdr("Alternative Hypothesis"))
  rightFrame <- tkframe(top)
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  pFrame <- tkframe(rightFrame)
  pVariable <- tclVar(".5")
  pField <- ttkentry(pFrame, width="6", textvariable=pVariable)
  radioButtons(name="test", buttons=c("normal", "corrected", "exact"), 
               labels=gettextRcmdr(c("Normal approximation", "Normal approximation with\ncontinuity correction", "Exact binomial")), 
               title=gettextRcmdr("Type of Test"))
  radioButtons(name="testVer", buttons=c("ordinaryVer", "wilsonVer"), 
               labels=gettextRcmdr(c("Ordinary (textbook default)", "Wilson score (R default)")), 
               title=gettextRcmdr("Approximation version"))
  
  tkgrid(labelRcmdr(pFrame, text=gettextRcmdr("Null hypothesis: p = "), fg="blue"), pField, sticky="nw")
  tkgrid(pFrame, sticky="nw")
  tkgrid(labelRcmdr(rightFrame, text=""))
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level: "), fg="blue"), confidenceField, sticky="nw")
  tkgrid(confidenceFrame, sticky="nw")
  tkgrid(alternativeFrame, rightFrame, sticky="nw")
  tkgrid(testFrame, testVerFrame, sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="nw")
  tkgrid.configure(confidenceField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}

##################################
# Partial least squares (and PCR)	
plsRegressionModel <- function(){
  initializeDialog(title=gettextRcmdr("Multivariate regression"))
  variablesFrame1 <- tkframe(top)
  variablesFrame2 <- tkframe(top)
  .numeric <- Numeric()
  .variable <- Variables()
  xBox <- variableListBox(variablesFrame1, .numeric, selectmode="multiple",
                          title=gettextRcmdr("Explanatory variables (pick one or more)"))
  yBox <- variableListBox(variablesFrame2, .variable, title=gettextRcmdr("Response variable (pick one)"))
  UpdateModelNumber()
  modelName <- tclVar(paste("MVRModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  compFrame <- tkframe(top)
  compVar <- tclVar("3")
  compEntry <- ttkentry(compFrame, width="3", textvariable=compVar)
  
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    .activeDataSet <- ActiveDataSet()
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=plsRegressionModel, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=plsRegressionModel, message=gettextRcmdr("No explanatory variables selected."))
      return()
    }
    if (is.element(y, x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=plsRegressionModel, message=gettextRcmdr("Response and explanatory variables must be different."))
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    ncomp <- tclvalue(compVar)
    if(trim.blanks(ncomp) == gettextRcmdr("")){
      ncomp <- 1
      warning('Number of components must be specified')
    }
    
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      UpdateModelNumber(-1)
      errorCondition(recall=plsRegressionModel, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        plsRegressionModel()
        return()
      }
    }
    
    validation <- as.character(tclvalue(validationVariable))
    if(validation == gettextRcmdr("LOO"))
      validate <- ", validation='LOO'"
    else {
      if(validation == gettextRcmdr("none"))
        validate <- ""
      else
        validate <- paste(", validation='CV', segments=", validation, sep="")
    }
    jackknife <- tclvalue(jackknifeVariable)
    if(jackknife == gettextRcmdr("1")){
      validate <- paste(validate, ", jackknife=TRUE", sep="")
      if(validation == gettextRcmdr("none")){
        UpdateModelNumber(-1)
        errorCondition(recall=plsRegressionModel, message=gettextRcmdr("Cannot perform jackknife without cross-validation."))
        return()
      }
    }
    
    pcrpls <- as.character(tclvalue(pcrplsrVariable))
    .the.y <- justDoIt(paste(".the.y <- ", ActiveDataSet(), "$", y, sep=""))
    if(is.factor(.the.y)){
      .the.new.y <- justDoIt(".the.new.y <- matrix(0, length(.the.y), length(levels(.the.y))-1)")
      for(i in 1:(length(levels(.the.y))-1)){
        .the.new.y[,i] <- justDoIt(paste(".the.new.y[,", i, "] <- (.the.y==levels(.the.y)[", i, "])*1",sep=""))
      }
      justDoIt(paste(ActiveDataSet(), "$", y, " <- .the.new.y", sep=""))
    }
    command <- paste(pcrpls, "(", y, "~", paste(x, collapse="+"),
                     ", data=", ActiveDataSet(), subset, ", ncomp=", as.numeric(ncomp), validate, ")", sep="")
    #    logger(paste(modelValue, " <- ", command, sep=""))
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    if(is.factor(.the.y)){
      justDoIt(paste(ActiveDataSet(), "$", y, " <- .the.y", sep=""))
      justDoIt("rm('.the.new.y')")
      doItAndPrint("#Dummy response")
    }
    justDoIt("rm('.the.y')")
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
    addComp <- tclvalue(addCompVariable)
    if (addComp == "1") {
      initializeDialog(subdialog, title=gettextRcmdr("Number of components"))
      tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Number of components to retain:"), fg="blue"), sticky="w")    
      sliderFrame <- tkframe(subdialog)
      sliderValue <- tclVar("1")
      componentsSlider <- tkscale(sliderFrame, from=1, to=ncomp, showvalue=FALSE, variable=sliderValue,
                                  resolution=1, orient="horizontal")
      componentsShow <- labelRcmdr(sliderFrame, textvariable=sliderValue, width=2, justify="right")
      onOKsub <- function() {
        closeDialog(subdialog)
        putRcmdr("ncomponents", as.numeric(tclvalue(sliderValue)))
      }
      subOKCancelHelp()
      tkgrid(componentsSlider, componentsShow, sticky="nw")
      tkgrid(sliderFrame, sticky="w")
      tkgrid(subButtonsFrame, sticky="w")
      dialogSuffix(subdialog, onOK=onOKsub, rows=2, columns=1, focus=subdialog)
      if ((ncomponents <- getRcmdr("ncomponents")) > 0){
        for(i in 1:ncomponents){
          var <- paste("Comp", i, sep="")
          if (is.element(var, Variables())) {
            if ("no" == tclvalue(checkReplace(var))) next
          }
          justDoIt(paste(.activeDataSet, "$Comp", i, " <- scores(", modelValue, ")[,", i, "]", sep=""))
          logger(paste(.activeDataSet, "$Comp", i, " <- scores(", modelValue, ")[,", i, "]", sep=""))
        }
        activeDataSet(.activeDataSet)
      }
    }
    activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="mvr", model=TRUE)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, row=1, column=1, columnspan=2, sticky="n")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame2, text="    "), sticky="nw")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame2, sticky="w", row=2, column=1)
  tkgrid(variablesFrame1, sticky="w", row=2, column=2)
  tkgrid(subsetFrame, row=3, column=1, sticky="w")
  tkgrid(labelRcmdr(compFrame, text=gettextRcmdr("Number of components")), compEntry, sticky="w")
  tkgrid(compFrame, row=3, column=2, sticky="w")
  radioButtons(name="validation", buttons=c("none", "LOO", "CV10", "CV5"), values=c("none", "LOO", "10", "5"),
               labels=gettextRcmdr(c("None", "Leave-one-out", "10-fold", "5-fold")), title=gettextRcmdr("Cross validation"))
  tkgrid(validationFrame, row=4, column=1, rowspan=2, columnspan=1, sticky="w")
  radioButtons(name="pcrplsr", buttons=c("pcr", "pls", "cppls"), values=c("pcr", "plsr", "cppls"), initialValue = "plsr",
               labels=gettextRcmdr(c("Principal components", "Partial least squares", "Canonical PLS")), title=gettextRcmdr("Type of regression"))
  tkgrid(pcrplsrFrame, row=4, column=2, rowspan=1, columnspan=1, sticky="w")
  checkBoxes(frame="optionsFrame", boxes=c("jackknife","addComp"), initialValues=c("0","0"),
             labels=gettextRcmdr(c("Jackknifing","Add scores to data set")))
  tkgrid(optionsFrame, row=5, column=2, columnspan=1, sticky="w")
  tkgrid(buttonsFrame, row=6, column=1, columnspan=2, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=6, columns=2)
}

########################
## Principal components
principalComponentsStat <- function(){
  initializeDialog(title=gettextRcmdr("Principal Components Analysis (JW)"))
  UpdateModelNumber()
  modelName <- tclVar(paste("PCAModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettextRcmdr("Variables (pick two or more)"))
  subsetBox()
  checkBoxes(frame="optionsFrame", boxes=c("center", "correlations", "screeplot", "addPC"), initialValues=c("1", "1", "0", "0"),
             labels=gettextRcmdr(c("Center predictors", "Analyze correlation matrix", "Screeplot", "Add principal components to data set")))
  onOK <- function(){ # Actions to perform
    putRcmdr("ncomponents", 0)
    x <- getSelection(xBox)
    nvar <- length(x)
    center <- tclvalue(centerVariable)
    correlations <- tclvalue(correlationsVariable)
    subset <- tclvalue(subsetVariable)
    screeplot <- tclvalue(screeplotVariable)
    addPC <- tclvalue(addPCVariable)
    closeDialog()
    if (2 > length(x)) {
      errorCondition(recall=principalComponentsStat, message=gettextRcmdr("Fewer than 2 variables selected."))
      return()
    }
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      UpdateModelNumber(-1)
      errorCondition(recall=principalComponentsStat, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        principalComponents()
        return()
      }
    }
    if (trim.blanks(subset) == "" || trim.blanks(subset) == gettextRcmdr("<all valid cases>")){
		subset <- ""
		putRcmdr("modelWithSubset", FALSE)
	} else {
		subset <- paste(", subset=", subset, sep="")
		putRcmdr("modelWithSubset", TRUE)
	}
    correlations <- if (correlations == "1") "TRUE" else "FALSE"
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("\n  ", modelValue, " <- prcomp(~", paste(x, collapse = "+"), 
                       ", scale.=", correlations, ", data=", .activeDataSet, 
                       subset, ")", sep = ""))
    cmds <- character(8)
    cmds[1] <- "local({"
    cmds[2] <- '  cat("\\nComponent loadings:\\n")'
    cmds[3] <- paste("  print(", modelValue, "$rotation)", sep="")
    cmds[4] <- '  cat("\\nComponent variances:\\n")'
    cmds[5] <- paste("  print(", modelValue, "$sdev^2)", sep="")
    cmds[6] <- '  cat("\\n")'
    cmds[7] <- paste("  print(summary(", modelValue, "))", sep="")
    cmds[8] <- if (screeplot == "1") paste("  screeplot(", modelValue, ")\n", sep="") else ""
    cmds <- paste(cmds[1:(7 + (screeplot == "1"))], collapse="\n")
    commands <- ""
    
    if (addPC == "1") {
      initializeDialog(subdialog, title=gettextRcmdr("Number of Components"))
      tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Number of components to retain:"), fg="blue"), sticky="w")    
      sliderFrame <- tkframe(subdialog)
      sliderValue <- tclVar("1")
      componentsSlider <- tkscale(sliderFrame, from=1, to=nvar, showvalue=FALSE, variable=sliderValue,
                                  resolution=1, orient="horizontal")
      componentsShow <- labelRcmdr(sliderFrame, textvariable=sliderValue, width=2, justify="right")
      onOKsub <- function() {
        closeDialog(subdialog)
        putRcmdr("ncomponents", as.numeric(tclvalue(sliderValue)))
      }
      subOKCancelHelp()
      tkgrid(componentsSlider, componentsShow, sticky="nw")
      tkgrid(sliderFrame, sticky="w")
      tkgrid(subButtonsFrame, sticky="w")
      dialogSuffix(subdialog, onOK=onOKsub, rows=2, columns=1, focus=subdialog, force.wait=TRUE)
      if ((ncomponents <- getRcmdr("ncomponents")) > 0){
        prefix <- paste("  ", .activeDataSet, " <<- within(", .activeDataSet, ", {", sep="")
        if (screeplot != "1") prefix <- paste("\n", prefix, sep="")
        commands <- character(ncomponents)
        for(i in 1:ncomponents){
          var <- paste("PC", i, sep="")
          if (is.element(var, Variables())) {
            if ("no" == tclvalue(checkReplace(var))) next
          }
          commands[ncomponents - i + 1] <- paste("    PC", i, " <- ", modelValue, "$x[,", 
                                                 i, "]", sep = "")
        }
        suffix <- "  })"
        commands <- paste(c(prefix, commands, suffix), collapse="\n")
      }
    }
    doItAndPrint(commands)
    doItAndPrint(paste(cmds, "\n})", sep=""))
    if (addPC == "1") activeDataSet(.activeDataSet, flushDialogMemory=FALSE)
    activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="prcomp")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="n")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=4, columns=1)
}

#########################################
# Linear/quadratic discriminant analysis
discriminantAnalysis <- function(){
  initializeDialog(title=gettextRcmdr("Discriminant analysis"))
  UpdateModelNumber()
  modelName <- tclVar(paste("DAModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  .factors <- Factors()
  xBox <- variableListBox(variablesFrame, .numeric, selectmode="multiple",
                          title=gettextRcmdr("Explanatory variables (pick one or more)"))
  yBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Response variable (pick one)"))
  cvFrame <- tkframe(top)
  priorFrame <- tkframe(top)
  predFrame <- tkframe(top)
  postFrame <- tkframe(top)
  subsetBox()
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=discriminantAnalysis, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=discriminantAnalysis, message=gettextRcmdr("No explanatory variables selected."))
      return()
    }
    if (is.element(y, x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=discriminantAnalysis, message=gettextRcmdr("Response and explanatory variables must be different."))
      return()
    }
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      UpdateModelNumber(-1)
      errorCondition(recall=discriminantAnalysis, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        discriminantAnalysis()
        return()
      }
    }
    the.subset <- tclvalue(subsetVariable)
    if (trim.blanks(the.subset) == gettextRcmdr("<all valid cases>") || trim.blanks(the.subset) == ""){
      subset <- ""
      the.subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", the.subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    linQuad <- as.character(tclvalue(linQuadVariable))
    if(linQuad == gettextRcmdr("linear")){
      command <- paste("lda(", y, "~", paste(x, collapse="+"),
                       ", data=", ActiveDataSet(), subset, "", sep="")
    }
    else{
      command <- paste("qda(", y, "~", paste(x, collapse="+"),
                       ", data=", ActiveDataSet(), subset, "", sep="")
    }
    prior <- as.character(tclvalue(priorVariable))
    if(prior == gettextRcmdr("equal")){
      try(eval(parse(text=paste("g <- length(levels(", ActiveDataSet(), "$", y, "))", sep=""))))
      command <- paste(command, ", prior=rep(", 1/g, ",", g, ")", sep="")
    }
    the.cv <- tclvalue(cvVariable)
    if(the.cv == gettextRcmdr("1")){
      command.cv <- paste(command, ", CV=TRUE)", sep="")
      command <- paste(command, ")", sep="")
    }
    else {
      command <- paste(command, ")", sep="")
    }
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    #    logger(paste(modelValue," <- ", command, sep=""))
    doItAndPrint(paste(modelValue," <- ", command, sep=""))
    if(the.cv == gettextRcmdr("0")){
      doItAndPrint(paste("confusion(", ActiveDataSet(), "$", y, "[", the.subset, "]", ", predict(", modelValue, ")$class)  # confusion matrix", sep=""))
    }
    else {
      doItAndPrint(paste("confusion(", ActiveDataSet(), "$", y, ", ", command.cv, "$class)  # confusion matrix", sep=""))
    }
    the.predDisp <- as.character(tclvalue(displayPredVariable))
    the.predSave <- as.character(tclvalue(savePredVariable))
    if(the.predDisp == gettextRcmdr("1")){
      if(the.cv == gettextRcmdr("0")){
        doItAndPrint(paste("predict(", modelValue, ")$class  # Predicted classes", sep=""))
      }
      else {
        doItAndPrint(paste(command.cv, "$class  # Predicted classes", sep=""))
      }
    }
    if(the.predSave == gettextRcmdr("1")){
      .activeDataSet <- ActiveDataSet()
      if(the.cv == gettextRcmdr("0")){
        justDoIt(paste(.activeDataSet, "$DA.class <- predict(", modelValue, ")$class", sep=""))
        logger(paste(.activeDataSet, "$DA.class <- predict(", modelValue, ")$class", sep=""))
      }
      else {
        justDoIt(paste(.activeDataSet, "$DA.class <- ", command.cv, "$class", sep=""))
        logger(paste(.activeDataSet, "$DA.class <- ", command.cv, "$class", sep=""))
      }
      activeDataSet(.activeDataSet)
    }
    the.postDisp <- as.character(tclvalue(displayPostVariable))
    the.postSave <- as.character(tclvalue(savePostVariable))
    if(the.postDisp == gettextRcmdr("1")){
      if(the.cv == gettextRcmdr("0")){
        doItAndPrint(paste("predict(", modelValue, ")$posterior  # Posterior probabilities", sep=""))
      }
      else {
        doItAndPrint(paste(command.cv, "$posterior  # Posterior probabilities", sep=""))
      }
    }
    if(the.postSave == gettextRcmdr("1")){
      .activeDataSet <- ActiveDataSet()
      if(the.cv == gettextRcmdr("0")){
        justDoIt(paste(.activeDataSet, "$DA.posterior <- predict(", modelValue, ")$posterior", sep=""))
        logger(paste(.activeDataSet, "$DA.posterior <- predict(", modelValue, ")$posterior", sep=""))
      }
      else {
        justDoIt(paste(.activeDataSet, "$DA.posterior <- ", command.cv, "$posterior", sep=""))
        logger(paste(.activeDataSet, "$DA.posterior <- ", command.cv, "$posterior", sep=""))
      }
      activeDataSet(.activeDataSet)
    }
    activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="lda")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, row=1, column=1, columnspan=2, sticky="n")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w", row=2, column=1, columnspan=2)
  tkgrid(subsetFrame, sticky="w", row=3, column=1, columnspan=1, rowspan=1)
  checkBoxes(frame="cvFrame", boxes=c("cv"), initialValues=c("0"),
             labels=gettextRcmdr(c("Cross-validation (LOO)")))
  tkgrid(cvFrame, row=3, column=2, columnspan=1, sticky="w")
  radioButtons(name="linQuad", buttons=c("Linear", "Quadratic"), values=c("linear", "quadratic"),
               labels=gettextRcmdr(c("Linear", "Quadratic")), title=gettextRcmdr("Type"))
  tkgrid(linQuadFrame, row=4, column=1, columnspan=1, sticky="w")
  radioButtons(name="prior", buttons=c("Empirical", "Equal"), values=c("emprirical", "equal"),
               labels=gettextRcmdr(c("Empirical", "Equal")), title=gettextRcmdr("Priors"))
  tkgrid(priorFrame, row=4, column=2, columnspan=1, sticky="w")
  checkBoxes(frame="predFrame", boxes=c("displayPred","savePred"), initialValues=c("0","0"),
             labels=gettextRcmdr(c("Display predictions","Save predictions")))
  tkgrid(predFrame, row=5, column=1, columnspan=1, sticky="w")
  checkBoxes(frame="postFrame", boxes=c("displayPost","savePost"), initialValues=c("0","0"),
             labels=gettextRcmdr(c("Display post. prob.","Save post. prob.")))
  tkgrid(postFrame, row=5, column=2, columnspan=1, sticky="w")
  tkgrid(buttonsFrame, row=6, column=1, columnspan=2, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=6, columns=2)
}

###################################################
# Hierarchical clustering of variables
hierarchicalClusterVariable <- function () {
	solutionNumber = length(listHclustSolutions())
	defaults <- list(initial.x = NULL, initial.clusMethod = "ward", initial.distance = "cor1",  
			initial.subset = gettextRcmdr ("<all valid cases>"),
			initial.dendro = 1, initial.tab=0)
	dialog.values <- getDialog("hierarchicalClusterVariable", defaults)
	initializeDialog(title = gettextRcmdr("Hierarchical Clustering of Variables"), use.tabs=TRUE)
	solutionFrame <- tkframe(dataTab)
	solutionName <- tclVar(paste("HClust.", (solutionNumber + 
								1), sep = ""))
	solutionField <- ttkentry(solutionFrame, width = "20", textvariable = solutionName)
	dataFrame <- tkframe(dataTab)
	xBox <- variableListBox(dataFrame, Numeric(), selectmode = "multiple", 
			title = gettextRcmdr("Variables (pick one or more)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	subsetBox(dataFrame, subset.expression = dialog.values$initial.subset)
	radioButtons(optionsTab, name = "method", buttons = c("ward", "single", 
					"complete", "average", "mcquitty", "median", "centroid"), 
			labels = gettextRcmdr(c("Ward's Method", "Single Linkage", 
							"Complete Linkage", "Average Linkage", "McQuitty's Method", 
							"Median Linkage", "Centroid Linkage")), title = gettextRcmdr("Clustering Method"), 
			initialValue = dialog.values$initial.clusMethod)
	optionsFrame <- tkframe(optionsTab)
	radioButtons(optionsFrame, name = "distanceType", buttons = c("cor1", "cor2"),
			labels = gettextRcmdr(c("Correlation", "Absolute correlation")), 
			title = gettextRcmdr("Distance Measure"), 
			initialValue = dialog.values$initial.distance)
	checkFrame <- tkframe(optionsFrame)
	plotDendro <- tclVar(dialog.values$initial.dendro)
	plotCB <- ttkcheckbutton(checkFrame)
	tkconfigure(plotCB, variable = plotDendro)
	onOK <- function() {
	    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
		x <- getSelection(xBox)
		nvar <- length(x)
		clusMethod <- tclvalue(methodVariable)
		distance <- tclvalue(distanceTypeVariable)
		subset <- trim.blanks(tclvalue(subsetVariable))
		dendro <- tclvalue(plotDendro)
		solution <- trim.blanks(tclvalue(solutionName))
		if (length(x) == 0) {
			errorCondition(recall = hierarchicalClusterVariable, message = gettextRcmdr("No variables selected."))
			return()
		}
		putDialog("hierarchicalClusterVariable", list(initial.x = x, initial.clusMethod = clusMethod, 
						initial.distance = distance, initial.subset = subset,
						initial.dendro = dendro, initial.tab = tab))
		closeDialog()
		varFormula <- paste(x, collapse = "+")
		vars <- paste(x, collapse = ",", sep = "")
		.activeDataSet <- ActiveDataSet()
		dset <- if (subset == gettextRcmdr("<all valid cases>")) 
					.activeDataSet
				else {
					paste(.activeDataSet, "[", .activeDataSet, "$", subset, 
							", ]", sep = "")
				}
		xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", 
				dset, ")", sep = "")
		if(distance=="cor1") {
		  dx <- paste("dist(1-cor(", xmat, "))", sep="")
		  distlab <- "correlation"
		}
		else {
		  dx <- paste("dist(1-abs(cor(", xmat, ")))", sep="")
		  distlab <- "absolute correlation"
		}
		command <- paste("hclust(", dx, " , method= ", "\"", 
				clusMethod, "\"", ")", sep = "")
        doItAndPrint(paste(solution, " <- ", command, sep = ""))
		if (dendro == "1") {
			justDoIt(paste("plot(", solution, ", main= ", "\"", 
							"Cluster Dendrogram for Solution ", solution, 
							"\"", ", xlab= ", "\"", "Observation Number in Data Set ", 
							dset, "\"", ", sub=", "\"", "Method=", clusMethod, 
							"; Distance=", distlab, "\"", ")", sep = ""))
			logger(paste("plot(", solution, ", main= ", "\"", 
							"Cluster Dendrogram for Solution ", solution, 
							"\"", ", xlab= ", "\"", "Observation Number in Data Set ", 
							dset, "\"", ", sub=", "\"", "Method=", clusMethod, 
							"; Distance=", distlab, "\"", ")", sep = ""))
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "hclust", reset = "hierarchicalClusterVariable", 
                 apply = "hierarchicalClusterVariable", model = TRUE)
	tkgrid(solutionField, sticky = "w")
	tkgrid(labelRcmdr(dataTab, text = gettextRcmdr("Clustering solution name:")), 
			solutionFrame, sticky = "w")
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(distanceTypeFrame, sticky = "w")
	tkgrid(labelRcmdr(checkFrame, text = "  "), sticky = "w")
	tkgrid(plotCB, labelRcmdr(checkFrame, text = gettextRcmdr("Plot Dendrogram  ")), 
			sticky = "w")
	tkgrid(checkFrame, sticky = "w")
    tkgrid(dataFrame, sticky="w")
    tkgrid(methodFrame, optionsFrame, sticky="nw")
	dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

############################
## Customized two-way table
enterTableNMBU <- function(){
  Library("vcd")
  env <- environment()
  initializeDialog(title=gettextRcmdr("Enter Two-Way Table"))
  outerTableFrame <- tkframe(top)
  assign(".tableFrame", tkframe(outerTableFrame), envir=env)
  setUpTable <- function(...){
    tkdestroy(get(".tableFrame", envir=env))
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    nrows <- as.numeric(tclvalue(rowsValue))
    ncols <- as.numeric(tclvalue(colsValue))
    make.col.names <- "labelRcmdr(.tableFrame, text='')"
    for (j in 1:ncols) {
      col.varname <- paste(".colname.", j, sep="")
      assign(col.varname, tclVar(j), envir=env)
      make.col.names <- paste(make.col.names, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                              col.varname, ")", sep="")
    }
    eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
    for (i in 1:nrows){
      varname <- paste(".tab.", i, ".1", sep="")
      assign(varname, tclVar("") , envir=env)
      row.varname <- paste(".rowname.", i, sep="")
      assign(row.varname, tclVar(i), envir=env)
      make.row <- paste("ttkentry(.tableFrame, width='5', textvariable=",
                        row.varname, ")", sep="")
      make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                        varname, ")", sep="")
      for (j in 2:ncols){
        varname <- paste(".tab.", i, ".", j, sep="")
        assign(varname, tclVar(""), envir=env)
        make.row <- paste(make.row, ", ", "ttkentry(.tableFrame, width='5', textvariable=",
                          varname, ")", sep="")
      }
      eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
    }
    tkgrid(get(".tableFrame", envir=env), sticky="w")
  }
  rowColFrame <- tkframe(top)
  rowsValue <- tclVar("2")
  rowsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=rowsValue,
                        resolution=1, orient="horizontal", command=setUpTable)
  rowsShow <- labelRcmdr(rowColFrame, textvariable=rowsValue, width=2, justify="right")
  colsValue <- tclVar("2")
  colsSlider <- tkscale(rowColFrame, from=2, to=10, showvalue=FALSE, variable=colsValue,
                        resolution=1, orient="horizontal", command=setUpTable)
  colsShow <- labelRcmdr(rowColFrame, textvariable=colsValue, width=2, justify="right")
  onOK <- function(){
    nrows <- as.numeric(tclvalue(rowsValue))
    ncols <- as.numeric(tclvalue(colsValue))
    cell <- 0
    counts <- rep(NA, nrows*ncols)
    row.names <- rep("", nrows)
    col.names <- rep("", ncols)
    for (i in 1:nrows) row.names[i] <-
      eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
    for (j in 1:ncols) col.names[j] <-
      eval(parse(text=paste("tclvalue(", paste(".colname.", j, sep=""),")", sep="")))
    for (i in 1:nrows){
      for (j in 1:ncols){
        cell <- cell+1
        varname <- paste(".tab.", i, ".", j, sep="")
        counts[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
      }
    }
    counts <- na.omit(counts)
    if (length(counts) != nrows*ncols){
      errorCondition(recall=enterTableNMBU, message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."), length(counts), nrows, ncols))
      return()
    }
    if (length(unique(row.names)) != nrows){
      errorCondition(recall=enterTableNMBU, message=gettextRcmdr("Row names are not unique."))
      return()
    }
    if (length(unique(col.names)) != ncols){
      errorCondition(recall=enterTableNMBU, message=gettextRcmdr("Column names are not unique."))
      return()
    }
    percents <- as.character(tclvalue(percentsVariable))
    chisq <- tclvalue(chisqVariable)
    chisqComp <- tclvalue(chisqComponentsVariable)
    expected <- tclvalue(expFreqVariable)
    fisher <- tclvalue(fisherVariable)
    closeDialog()
    command <- paste("matrix(c(", paste(counts, collapse=","), "), ", nrows, ", ", ncols,
                     ", byrow=TRUE)", sep="")
    
    commandA <- paste(command, "\n  dimnames(.Table) <- list(c(",paste(paste("'", row.names, "'", sep=""), collapse=", "), "), c(",paste(paste("'", col.names, "'", sep=""), collapse=", "), "))", sep="")
    
    command <- paste("local({\n  .Table <- ", commandA, '\n  cat("\\nFrequency table:\\n")\n  print(.Table)', sep="")
    command.2 <- paste("local({\n  .warn <- options(warn=-1)\n  .Table <- ", commandA, sep="")
    if (percents == "row") 
      command <- paste(command, '\n  cat("\\nRow percentages:\\n")\n  print(rowPercents(.Table))',
                       sep="")
    else if (percents == "column") 
      command <-  paste(command, '\n  cat("\\nColumn percentages:\\n")\n  print(colPercents(.Table))',
                        sep="")
    else if (percents == "total") 
      command <- paste(command, '\n  cat("\\nTotal percentages:\\n")\n  print(totPercents(.Table))',
                       sep="")
    if (chisq == 1) {
      command <- paste(command, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
      command.2 <- paste(command.2, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
      command <- paste(command, "\n  print(assocstats(.Table))", sep="")
      if (expected == 1)
        command <- paste(command, '\n  cat("\\nExpected counts:\\n")\n  print(.Test$expected)', sep="")
      if (chisqComp == 1) {
        command <- paste(command, '\n  cat("\\nChi-square components:\\n")\n  print(round(.Test$residuals^2, 2))', sep="")
        command <- paste(command, '\n  cat("\\nAdjusted residuals:\\n\")\n  ', "print(round((.Table-.Test$expected)/sqrt(.Test$expected*tcrossprod((1-apply(.Table,1,sum)/sum(.Table)),(1-apply(.Table,2,sum)/sum(.Table)))),2))", sep="")
      }
    }
    if (fisher == 1) command <- paste(command, "\n  print(fisher.test(.Table))")
    command <- paste(command, ' \n})', sep="")
    doItAndPrint(command)
    if (chisq == 1){
      command.2 <- paste(command.2, "\nputRcmdr('.expected.counts', .Test$expected)\n  options(.warn)\n})")
      justDoIt(command.2)
      warnText <- NULL
      expected <- getRcmdr(".expected.counts")
      if (0 < (nlt1 <- sum(expected < 1))) warnText <- paste(nlt1,
                                                             gettextRcmdr("expected frequencies are less than 1"))
      if (0 < (nlt5 <- sum(expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                                                             gettextRcmdr(" expected frequencies are less than 5"), sep="")
      if (!is.null(warnText)) Message(message=warnText,
                                      type="warning")
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="chisq.test")
  radioButtons(name="percents", buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"), values=c("row", "column", "total", "none"),
               initialValue="none", labels=gettextRcmdr(c("Row percentages", "Column percentages",  "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
  checkBoxes(frame="testsFrame", boxes=c("chisq", "chisqComponents", "expFreq", "fisher"), initialValues=c("1", "1", "1", "0"),
             labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
                                   "Print expected frequencies", "Fisher's exact test")))
  tkgrid(labelRcmdr(rowColFrame, text=gettextRcmdr("Number of Rows:")), rowsSlider, rowsShow, sticky="w")
  tkgrid(labelRcmdr(rowColFrame, text=gettextRcmdr("Number of Columns:")), colsSlider, colsShow, sticky="w")
  tkgrid(rowColFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter counts:"), fg="blue"), sticky="w")
  tkgrid(outerTableFrame, sticky="w")
  tkgrid(percentsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
  tkgrid(testsFrame, sticky="w")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=7, columns=2)
}

#############################
## Customized two-way table
twoWayTableNMBU <- function(){
  Library("vcd")
  initializeDialog(title=gettextRcmdr("Two-Way Table"))
  variablesFrame <- tkframe(top)
  .factors <- Factors()
  .numeric <- Numeric()
  countBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Counts (pick one)"))
  rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"))
  columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"))
  onOK <- function(){
    count <- getSelection(countBox)
    row <- getSelection(rowBox)
    column <- getSelection(columnBox)
    if (length(count) == 0 || length(row) == 0 || length(column) == 0){
      errorCondition(recall=twoWayTableNMBU, message=gettextRcmdr("You must select three variables."))
      return()
    }
    if (row == column) {
      errorCondition(recall=twoWayTableNMBU, message=gettextRcmdr("Row and column variables are the same."))
      return()
    }
    percents <- as.character(tclvalue(percentsVariable))
    chisq <- tclvalue(chisqTestVariable)
    chisqComp <- tclvalue(chisqComponentsVariable)
    expected <- tclvalue(expFreqVariable)
    fisher <- tclvalue(fisherTestVariable)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    commandA <- paste("table(data.frame(", row, "=rep(", .activeDataSet, "[,'", row, "'],", .activeDataSet, "[,'", count, "']),", column, "=rep(", .activeDataSet, "[,'", column, "'],", .activeDataSet, "[,'", count, "'])))", sep="")
    
    command <- paste("local({\n  .Table <- ", commandA, '\n  cat("\\nFrequency table:\\n")\n  print(.Table)', sep="")
    command.2 <- paste("local({\n  .warn <- options(warn=-1)\n  .Table <- ", commandA, sep="")
    if (percents == "row") 
      command <- paste(command, '\n  cat("\\nRow percentages:\\n")\n  print(rowPercents(.Table))',
                       sep="")
    else if (percents == "column") 
      command <-  paste(command, '\n  cat("\\nColumn percentages:\\n")\n  print(colPercents(.Table))',
                        sep="")
    else if (percents == "total") 
      command <- paste(command, '\n  cat("\\nTotal percentages:\\n")\n  print(totPercents(.Table))',
                       sep="")
    if (chisq == 1) {
      command <- paste(command, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
      command.2 <- paste(command.2, "\n  .Test <- chisq.test(.Table, correct=FALSE)", sep="")
      command <- paste(command, "\n  print(assocstats(.Table))", sep="")
      if (expected == 1)
        command <- paste(command, '\n  cat("\\nExpected counts:\\n")\n  print(.Test$expected)', sep="")
      if (chisqComp == 1) {
        command <- paste(command, '\n  cat("\\nChi-square components:\\n")\n  print(round(.Test$residuals^2, 2))', sep="")
        command <- paste(command, '\n  cat("\\nAdjusted residuals:\\n\")\n  ', "print(round((.Table-.Test$expected)/sqrt(.Test$expected*tcrossprod((1-apply(.Table,1,sum)/sum(.Table)),(1-apply(.Table,2,sum)/sum(.Table)))),2))", sep="")
      }
    }
    if (fisher == 1) command <- paste(command, "\n  print(fisher.test(.Table))")
    command <- paste(command, ' \n})', sep="")
    doItAndPrint(command)
    if (chisq == 1){
      command.2 <- paste(command.2, "\nputRcmdr('.expected.counts', .Test$expected)\n  options(.warn)\n})")
      justDoIt(command.2)
      warnText <- NULL
      expected <- getRcmdr(".expected.counts")
      if (0 < (nlt1 <- sum(expected < 1))) warnText <- paste(nlt1,
                                                             gettextRcmdr("expected frequencies are less than 1"))
      if (0 < (nlt5 <- sum(expected < 5))) warnText <- paste(warnText, "\n", nlt5,
                                                             gettextRcmdr(" expected frequencies are less than 5"), sep="")
      if (!is.null(warnText)) Message(message=warnText,
                                      type="warning")
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="xtabs")
  radioButtons(name="percents",
               buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
               values=c("row", "column", "total", "none"), initialValue="none",
               labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
  checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("1", "1", "1", "0"),
             labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
                                   "Print expected frequencies", "Fisher's exact test")))
  tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), labelRcmdr(variablesFrame, text="    "), getFrame(countBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(percentsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
  tkgrid(testsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=6, columns=1)
}

################################
# Analysis of variance
linearModelNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Linear model (regression/ANOVA) - specify model"))
  defaults <- list(initial.contr="contr.sum", initial.output=c("1","0","0","0","0"), initial.subset = gettextRcmdr("<all valid cases>"), initial.weight = gettextRcmdr("<no variable selected>"))
  dialog.values <- getDialog("linearModelNMBU", defaults)
  .activeModel <- ActiveModel()
  variables <- Variables()
  factors <- Factors()
  # To be able to recreate settings from former models many things are defined here that have to do with randomness and such
  chosen.factors <- FALSE
  .factorsLabel <- tclVar("Set factors")
  onFactors <- function(){
    initializeDialog(subdialog,title=gettextRcmdr("Set factors"))
    .variable <- Variables()
    yBox <- variableListBox(subdialog, .variable, title=gettextRcmdr("Convert to factor(s) (pick zero or more)"), selectmode="multiple")
    onOKsub <- function(){
      chosen <- getSelection(yBox)
      if (length(chosen) == 0){
        assign("chosen.factors", FALSE, envir=env)
        tclvalue(.factorsLabel) <- "Set factors"
        tkconfigure(factorsButton, foreground="black")
        if (GrabFocus()) tkgrab.release(subdialog)
        tkdestroy(subdialog)
        tkwm.deiconify(top)
        if (GrabFocus()) tkgrab.set(top)
        tkfocus(top)
        tkwait.window(top)
        return()
      }
      assign("chosen.factors", chosen, envir=env)
      tclvalue(.factorsLabel) <- "Factors set"
      tkconfigure(factorsButton, foreground="blue")
      tkdestroy(subdialog)
      tkwm.deiconify(top)
      if (GrabFocus()) tkgrab.set(top)
      tkfocus(top)
      tkwait.window(top)
    }
    subOKCancelHelp()
    tkgrid(getFrame(yBox), sticky="nw")
    tkgrid(subButtonsFrame, columnspan=2, sticky="w")
    dialogSuffix(subdialog, rows=2, columns=2, onOK=onOKsub, focus=subdialog)
  }
  
  groupsFrame <- tkframe(top)
  factorsButton <- tkbutton(groupsFrame, textvariable=.factorsLabel, command=onFactors, borderwidth=3)
  checkBoxes(frame="optionsFrame", boxes=c("reg","type1","type2","type3","typeR"), initialValues=dialog.values$initial.output, #c("1","0","0","0","0"),
             labels=gettextRcmdr(c("Regression","ANOVA 'type I test' (sequential)", "ANOVA 'type II test' (obeying marginality)", "ANOVA 'type III test' (ignoring marginality)", "ANOVA for regression")))
  onRandom <- function() {
    if (GrabFocus() && .Platform$OS.type != "windows") tkgrab.release(window)
    if (as.numeric(R.Version()$major) >= 2) print(help("mixlm"))
    else help("mixlm")
  }
  helpButton <- buttonRcmdr(optionsFrame, text=gettextRcmdr("Random effects help"), width="20", command=onRandom, borderwidth=3,
        image="::image::helpIcon", compound="left")
  tkgrid(helpButton, sticky="w")
  
  currentModel <- if (!is.null(.activeModel))
    class(get(.activeModel, envir=.GlobalEnv))[1] == "lm" || class(get(.activeModel, envir=.GlobalEnv))[1] == "lmm" || class(get(.activeModel, envir=.GlobalEnv))[1] == "mer"
  else FALSE
  if (currentModel) {
    currentFields <- formulaFields2(get(.activeModel, envir=.GlobalEnv))
    if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
    the.class <- class(get(.activeModel, envir=.GlobalEnv))
    if(grepl("~",currentFields$lhs,fixed=TRUE)){
      tmp <- strsplit(currentFields$lhs, "~", fixed=TRUE)[[1]]
      currentFields$lhs <- tmp[1]; currentFields$rhs <- tmp[2]
    }
  }
  if (isTRUE(getRcmdr("reset.model"))) {
    currentModel <- FALSE
    putRcmdr("reset.model", FALSE)
  }
  UpdateModelNumber()
  modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  onOK <- function(){
    modelValue <- trim.blanks(tclvalue(modelName))
    closeDialog()
    if (!is.valid.name(modelValue)){
      errorCondition(recall=linearModelNMBU, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
    subset <- tclvalue(subsetVariable)
    weight.var <- getSelection(weightComboBox)
    putDialog ("linearModelNMBU", list (initial.contr = tclvalue(contrVariable), 	
			   initial.output = c(tclvalue(regVariable),
               tclvalue(type1Variable), tclvalue(type2Variable), tclvalue(type3Variable), tclvalue(typeRVariable)), initial.subset = subset, initial.weight = weight.var))
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    weights <- if (weight.var == gettextRcmdr("<no variable selected>")) ""
    else paste(", weights=", weight.var, sep="")
    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=linearModelNMBU, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
      return()
    }
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=linearModelNMBU, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
      return()
    }
    if (is.element(modelValue, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        linearModel()
        return()
      }
    }
	old.contr <- options("contrasts")
    if(tclvalue(contrVariable) != gettextRcmdr("reml")){
      if(tclvalue(contrVariable)==gettextRcmdr("contr.sum")){
        command <- "options(contrasts=c('contr.sum','contr.poly'))"
      } else {
        command <- "options(contrasts=c('contr.treatment','contr.poly'))"
      }
      if(options("contrasts")$contrasts[1] != tclvalue(contrVariable))
        doItAndPrint(command)
    }
    .activeDataSet <- ActiveDataSet()
    if(!is.logical(chosen.factors)){
      variables <- Variables()
      command <- paste(.activeDataSet, ".tmp <- data.frame(", paste(variables[!is.element(variables,chosen.factors)], "=", ActiveDataSet(), "$", variables[!is.element(variables,chosen.factors)], ", ", sep="", collapse=""), paste(variables[is.element(variables,chosen.factors)], "=factor(", ActiveDataSet(), "$", variables[is.element(variables,chosen.factors)], "), ", sep="", collapse=""), sep="")
      command <- paste(substr(command,1,nchar(command)-2),")", sep="")
      doItAndPrint(command)
      .activeDataSet <- paste(.activeDataSet, ".tmp", sep="")
      activeDataSet(.activeDataSet)
    }
    
    type <- as.character(tclvalue(contrVariable))
	formula1 <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
    
    formula1 <- paste("lm(", formula1, sep="")
    command <- paste(formula1, ", data=", .activeDataSet, subset, weights, sep="")
    if(tclvalue(contrVariable)==gettextRcmdr("reml"))
      command <- paste(command, ", REML=TRUE)", sep="")
    else
      command <- paste(command, ")", sep="")
    #	logger(paste(modelValue, " <- ", command, sep=""))
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
	options(contrasts = old.contr$contrasts)	
    activeModel(modelValue)
	if(!is.null(ActiveModel())){
		if(tclvalue(regVariable)   == gettextRcmdr("1")) doItAndPrint(paste("summary(", modelValue, ")", sep=""))
		if(tclvalue(type1Variable) == gettextRcmdr("1")) doItAndPrint(paste("anova(", ActiveModel(), ")", sep=""))
		if(tclvalue(type2Variable) == gettextRcmdr("1")) doItAndPrint(paste("Anova(", ActiveModel(), ', type="II")', sep=""))
		if(tclvalue(type3Variable) == gettextRcmdr("1")) doItAndPrint(paste("Anova(", ActiveModel(), ', type="III")', sep=""))
		if(tclvalue(typeRVariable) == gettextRcmdr("1")) doItAndPrint(paste("anova_reg(", ActiveModel(), ')', sep=""))
	}
    #    if(tclvalue(mixVariable)   == gettextRcmdr("1")) mixed.modelGUI()
    tkfocus(CommanderWindow())
  }
  env <- environment()
  tkgrid(labelRcmdr(groupsFrame, text="    "), factorsButton, sticky="w")
  
  OKCancelHelp(helpSubject="linearModel", model=TRUE, reset="resetLinearModelNMBU", apply="linearModelNMBU")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", column=1, row=1, columnspan=2)
  modelFormula3()#.variables=variables, .factors=factors)
  
  subsetWeightFrame <- tkframe(top)

  subsetBox(window=subsetWeightFrame, model=TRUE)#subset.expression = dialog.values$initial.subset)
  weightComboBox <- variableComboBox(subsetWeightFrame, variableList=Numeric(), 
                                     initialSelection=dialog.values$initial.weight,
                                     title=gettextRcmdr("Weights"))
  tkgrid(getFrame(xBox), sticky="w", column=1, row=2, columnspan=1)
  tkgrid(groupsFrame, sticky="w", column=2, row=2, columnspan=1)
  tkgrid(outerOperatorsFrame, sticky="w", column=1, row=3, columnspan=2)
  tkgrid(formulaFrame, sticky="w", column=1, row=4, columnspan=2)
  tkgrid(subsetFrame, tklabel(subsetWeightFrame, text="   "),
         getFrame(weightComboBox), sticky="nw")
  radioButtons(name="contr", buttons=c("contr.sum", "contr.treatment", "reml"), values=c("contr.sum", "contr.treatment", "reml"), initialValue = dialog.values$initial.contr, #"contr.sum",
               labels=gettextRcmdr(c("Sum to zero (contr.sum)", "Reference level (contr.treatment)", "REML")), title=gettextRcmdr("Parameterization"))
  tkgrid(contrFrame, row=5, column=1, rowspan=1, columnspan=1, sticky="w")
  tkgrid(subsetWeightFrame, sticky="w", column=1, row=6)
  tkgrid(optionsFrame, row=5, column=2, columnspan=1, rowspan=2, sticky="w")
  tkgrid(buttonsFrame, sticky="w", column=1, row=7, columnspan=2)
  dialogSuffix(rows=7, columns=2,focus=lhsEntry, preventDoubleClick=TRUE)
}
resetLinearModelNMBU <- function(){
  putRcmdr("reset.model", TRUE)
  putDialog("linearModelNMBU", NULL)
  linearModelNMBU()
}


################################
# Customized GLM
generalizedLinearModelNMBU <- function(){
  families <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                "quasibinomial", "quasipoisson")
  links <- c("identity", "inverse", "log", "logit", "probit",
             "cloglog", "sqrt", "1/mu^2")
  availableLinks <- matrix(c(
    TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
    TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE,
    TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
    TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE,
    FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
    TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE),
    7, 8, byrow=TRUE)
  rownames(availableLinks) <- families
  colnames(availableLinks) <- links
  canonicalLinks <- c("identity", "logit", "log", "inverse", "1/mu^2", "logit", "log")
  names(canonicalLinks) <- families
  defaults <- list(initial.weight = gettextRcmdr("<no variable selected>"),initial.offset = gettextRcmdr("<no variable selected>"))
  dialog.values <- getDialog("generalizedLinearModelNMBU", defaults)
  initializeDialog(title=gettextRcmdr("Generalized Linear Model"))
  .activeModel <- ActiveModel()
  variables <- Variables()
  factors <- Factors()
  # To be able to recreate settings from former models many things are defined here that have to do with randomness and such
  chosen.factors <- FALSE
  .factorsLabel <- tclVar("Set factors")
  onFactors <- function(){
    initializeDialog(subdialog,title=gettextRcmdr("Set factors"))
    .variable <- Variables()
    yBox <- variableListBox(subdialog, .variable, title=gettextRcmdr("Convert to factor(s) (pick zero or more)"), selectmode="multiple")
    onOKsub <- function(){
      chosen <- getSelection(yBox)
      if (length(chosen) == 0){
        assign("chosen.factors", FALSE, envir=env)
        tclvalue(.factorsLabel) <- "Set factors"
        tkconfigure(factorsButton, foreground="black")
        if (GrabFocus()) tkgrab.release(subdialog)
        tkdestroy(subdialog)
        tkwm.deiconify(top)
        if (GrabFocus()) tkgrab.set(top)
        tkfocus(top)
        tkwait.window(top)
        return()
      }
      assign("chosen.factors", chosen, envir=env)
      tclvalue(.factorsLabel) <- "Factors set"
      tkconfigure(factorsButton, foreground="blue")
      tkdestroy(subdialog)
      tkwm.deiconify(top)
      if (GrabFocus()) tkgrab.set(top)
      tkfocus(top)
      tkwait.window(top)
    }
    subOKCancelHelp()
    tkgrid(getFrame(yBox), sticky="nw")
    tkgrid(subButtonsFrame, columnspan=2, sticky="w")
    dialogSuffix(subdialog, rows=2, columns=2, onOK=onOKsub, focus=subdialog)
  }
  
  groupsFrame <- tkframe(top)
  factorsButton <- tkbutton(groupsFrame, textvariable=.factorsLabel, command=onFactors, borderwidth=3)
  
  currentModel <- if (!is.null(.activeModel))
    class(get(.activeModel, envir=.GlobalEnv))[1] == "glm"
  else FALSE
  if (currentModel) {
    currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), glm=TRUE)
    if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
  }
  if (isTRUE(getRcmdr("reset.model"))) {
    currentModel <- FALSE
    putRcmdr("reset.model", FALSE)
  }
  modelFormula3()#.variables=variables, .factors=factors)
  UpdateModelNumber()
  modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  linkFamilyFrame <- tkframe(top)
  familyFrame <- tkframe(linkFamilyFrame)
  max.height <- getRcmdr("variable.list.height")
  familyBox <- tklistbox(familyFrame, height=min(max.height, length(families)),
						 exportselection="FALSE",
                         selectmode="single", background="white")
  familyScroll <- ttkscrollbar(familyFrame,
                               command=function(...) tkyview(familyBox, ...))
  tkconfigure(familyBox, yscrollcommand=function(...) tkset(familyScroll, ...))
  for (fam in families) tkinsert(familyBox, "end", fam)
  linkFrame <- tkframe(linkFamilyFrame)
  linkBox <- tklistbox(linkFrame, height=max.height, exportselection="FALSE",
                       selectmode="single", background="white")
 # subsetFrame <- tkframe(top)
  weightFrame <- tkframe(top)
  offsetFrame <- tkframe(top)
  subsetBox(model=TRUE)
  weightComboBox <- variableComboBox(weightFrame, variableList=Numeric(), 
                                     initialSelection=dialog.values$initial.weight,
                                     title=gettextRcmdr("Weights"))
  offsetComboBox <- variableComboBox(offsetFrame, variableList=Numeric(), 
                                     initialSelection=dialog.values$initial.offset,
                                     title=gettextRcmdr("Offset"))
  onFamilySelect <- function(){
    family <- families[as.numeric(tkcurselection(familyBox)) + 1]
    availLinks <- links[availableLinks[family,]]
    tkdelete(linkBox, "0", "end")
    for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
    canLink <- canonicalLinks[family]
    tkconfigure(linkBox, height=length(availLinks))
    tkselection.set(linkBox, which(canLink == availLinks) - 1)
  }
  onOK <- function(){
    .activeDataSet <- ActiveDataSet()
    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=generalizedLinearModelNMBU, model=TRUE, message=gettextRcmdr("Left-hand side of model empty."))
      return()
    }
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=generalizedLinearModelNMBU, model=TRUE, message=gettextRcmdr("Right-hand side of model empty."))
      return()
    }
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      errorCondition(recall=generalizedLinearModelNMBU, model=TRUE, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listGeneralizedLinearModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        closeDialog()
        generalizedLinearModelNMBU()
        return()
      }
    }
#	old.contr <- options("contrasts")
#    if(tclvalue(contrVariable)==gettextRcmdr("contr.sum")){
#      command <- "options(contrasts=c('contr.sum','contr.poly'))"
#    } else {
#      command <- "options(contrasts=c('contr.treatment','contr.poly'))"
#    }
#    if(options("contrasts")$contrasts[1] != tclvalue(contrVariable))
#      doItAndPrint(command)
    if(!is.logical(chosen.factors)){
      variables <- Variables()
      command <- paste(.activeDataSet, ".tmp <- data.frame(", paste(variables[!is.element(variables,chosen.factors)], "=", ActiveDataSet(), "$", variables[!is.element(variables,chosen.factors)], ", ", sep="", collapse=""), paste(variables[is.element(variables,chosen.factors)], "=factor(", ActiveDataSet(), "$", variables[is.element(variables,chosen.factors)], "), ", sep="", collapse=""), sep="")
      command <- paste(substr(command,1,nchar(command)-2),")", sep="")
      doItAndPrint(command)
      .activeDataSet <- paste(.activeDataSet, ".tmp", sep="")
      activeDataSet(.activeDataSet)
    }
    formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
    family <- families[as.numeric(tkcurselection(familyBox)) + 1]
    availLinks <- links[availableLinks[family,]]
    link <- availLinks[as.numeric(tkcurselection(linkBox)) + 1]
    subset <- tclvalue(subsetVariable)
    closeDialog()
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    weight.var <- getSelection(weightComboBox)
    offset.var <- getSelection(offsetComboBox)
    putDialog ("generalizedLinearModelNMBU", list (initial.offset = offset.var,initial.weight = weight.var))
    weights <- if (weight.var == gettextRcmdr("<no variable selected>")) ""
    else paste(", weights=", weight.var, sep="")
    offset <- if (offset.var == gettextRcmdr("<no variable selected>")) ""
    else paste(", offset=", offset.var, sep="")
    command <- paste("glm(", formula, ", family=", family, "(", link,
                     "), data=", ActiveDataSet(), subset, weights, offset, ")", sep="")
    #    logger(paste(modelValue, " <- ", command, sep=""))
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
    activeModel(modelValue)
    model.class <- eval(parse(text=paste("class(",modelValue,")",sep="")))
    if(model.class[1]=="glmerMod" || model.class[1]=="lmerMod"){
      doItAndPrint(paste("Anova(", modelValue, ",test='Chisq',type=3)", sep=""))
    } else {
      doItAndPrint(paste("Anova(", modelValue, ",test='LR',type=3)", sep=""))
    }
    doItAndPrint(paste("logLik(", modelValue, ")",sep=""))
#	options(contrasts = old.contr$contrasts)
    tkfocus(CommanderWindow())
  }
  env <- environment()
  tkgrid(labelRcmdr(groupsFrame, text="    "), factorsButton, sticky="w")
  
  OKCancelHelp(helpSubject="generalizedLinearModel", model=TRUE, reset="resetGLMNMBU", apply="generalizedLinearModelNMBU")
  helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", row=1, column=1, columnspan=3)
  tkgrid(getFrame(xBox), sticky="w", row=2, column=1, columnspan=1)
  tkgrid(groupsFrame, sticky="w", column=2, row=2, columnspan=1)
  tkgrid(outerOperatorsFrame, sticky="w", row=3, column=1, columnspan=3)
  tkgrid(formulaFrame, sticky="w", row=4, column=1, columnspan=3)
#  tkgrid(labelRcmdr(subsetFrame, text=gettextRcmdr("Subset expression")), getFrame(subsetBox), sticky="w")
  tkgrid(labelRcmdr(weightFrame, text=gettextRcmdr("")), getFrame(weightComboBox), sticky="w")
  tkgrid(labelRcmdr(offsetFrame, text=gettextRcmdr("")), getFrame(offsetComboBox), sticky="w")
  tkgrid(subsetFrame, sticky="w", row=5, column=1, columnspan=1)
  tkgrid(weightFrame, sticky="w", row=5, column=2, columnspan=1)
  tkgrid(offsetFrame, sticky="w", row=5, column=3, columnspan=1)
#  radioButtons(name="contr", buttons=c("contr.sum", "contr.treatment"), values=c("contr.sum", "contr.treatment"), initialValue = dialog.values$initial.contr, #"contr.sum",
#               labels=gettextRcmdr(c("Sum to zero (contr.sum)", "Reference level (contr.treatment)")), title=gettextRcmdr("Parameterization"))
#  tkgrid(contrFrame, row=5, column=1, rowspan=1, columnspan=1, sticky="w")
#  spaceFrame <- tkframe(top)
#  tkgrid(labelRcmdr(spaceFrame, text=gettextRcmdr(" ")), sticky="w")
#  tkgrid(spaceFrame, sticky="w", row=6, column=1, columnspan=2)
#  tkgrid(labelRcmdr(weightsFrame, text=gettextRcmdr("Weights (optional):")), weights, sticky="w")
#  tkgrid(weightsFrame, sticky="w", row=7, column=1, columnspan=1)
#  tkgrid(labelRcmdr(offsetFrame, text=gettextRcmdr("    Offset (optional):")), offset, sticky="w")
#  tkgrid(offsetFrame, sticky="w", row=7, column=2, columnspan=1)
#  tkgrid(subsetFrame, sticky="w", row=8, column=1, columnspan=2)
  tkgrid(labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Family (double-click to select)"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"),
         labelRcmdr(linkFamilyFrame, text="   "), labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Link function"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
  tkgrid(familyBox, familyScroll, sticky="nw")
  tkgrid(linkBox, sticky="nw")
  tkgrid(familyFrame, labelRcmdr(linkFamilyFrame, text="   "), linkFrame, sticky="nw")
  tkgrid(linkFamilyFrame, sticky="w", row=9, column=1, columnspan=3)
  tkgrid(buttonsFrame, sticky="w", row=10, column=1, columnspan=3)
  tkgrid.configure(familyScroll, sticky="ns")
  fam <- if (currentModel) which(currentFields$family == families) - 1
  else 1
  tkselection.set(familyBox, fam)
  availLinks <- links[availableLinks[fam + 1,]]
  for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
  tkconfigure(linkBox, height=length(availLinks))
  lnk <- if (currentModel) which(currentFields$link == availLinks) - 1
  else 0
  tkselection.set(linkBox, lnk)
  tkbind(familyBox, "<Double-ButtonPress-1>", onFamilySelect)
  dialogSuffix(rows=10, columns=3, focus=lhsEntry, preventDoubleClick=TRUE)
}
resetGLMNMBU <- function(){
  putRcmdr("reset.model", TRUE)
  putDialog("generalizedLinearModelNMBU", NULL)
  putDialog("generalizedLinearModelNMBU", NULL, resettable=FALSE)
  generalizedLinearModelNMBU()
}


################################
# Customized multiLogit
multinomialLogitModelNMBU <- function(){
  Library("nnet")
  initializeDialog(title=gettextRcmdr("Multinomial Logit Model"))
  .activeModel <- ActiveModel()
  .activeDataSet <- ActiveDataSet()
  currentModel <- if (!is.null(.activeModel))
    class(get(.activeModel, envir=.GlobalEnv))[1] == "multinom"
		else FALSE
  if (currentModel) {
    currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
    if (currentFields$data != .activeDataSet) currentModel <- FALSE
  }
	if (isTRUE(getRcmdr("reset.model"))) {
		currentModel <- FALSE
		putRcmdr("reset.model", FALSE)
	}
  UpdateModelNumber()
  modelName <- tclVar(paste("MLM.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  baseLevel <- tclVar("")
  baseFrame <- tkframe(top)
  baseV <- ttkentry(baseFrame, width="20", textvariable=baseLevel)
#  weightsName <- tclVar(gettextRcmdr("<none>"))
  weightsFrame <- tkframe(top)
  weightComboBox <- variableComboBox(weightsFrame, variableList=Numeric(), 
                                     initialSelection=gettextRcmdr("<no variable selected>"),
                                     title=gettextRcmdr("Weights"))
#  weights <- ttkentry(weightsFrame, width="20", textvariable=weightsName)
  
  onOK <- function(){
    modelValue <- trim.blanks(tclvalue(modelName))
    baseValue <- trim.blanks(tclvalue(baseLevel))
#    weightsValue <- tclvalue(weightsName)
    weight.var <- getSelection(weightComboBox)
    weights <- if (weight.var == gettextRcmdr("<no variable selected>")) ""
    else paste(", weights=", weight.var, sep="")
    closeDialog()
    if (!is.valid.name(modelValue)){
      errorCondition(recall=multinomialLogitModelNMBU, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=multinomialLogitModelNMBU, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
      return()
    }
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=multinomialLogitModelNMBU, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
      return()
    }
    if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=get(.activeDataSet, envir=.GlobalEnv)))){
      #        if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=eval(parse(text=.activeDataSet), envir=.GlobalEnv)))){
      errorCondition(recall=multinomialLogitModelNMBU, message=gettextRcmdr("Response variable must be a factor"))
      return()
    }
    if (baseValue!="" && !is.element(baseValue, justDoIt(paste("levels(", ActiveDataSet(), "$", tclvalue(lhsVariable), ")", sep="")))){
      errorCondition(recall=multinomialLogitModelNMBU, message=gettextRcmdr("'Base level' must be a level used in the response."))
      return()
    }
    if (is.element(modelValue, listMultinomialLogitModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        multinomialLogitModelNMBU()
        return()
      }
    }
    if(baseValue != ""){
      command <- paste("levels(", .activeDataSet, "$", tclvalue(lhsVariable), ")", sep="")
      effLevsO <- justDoIt(command)
      effLevs <- c(baseValue, effLevsO[!is.element(effLevsO,baseValue)])
      command <- paste(.activeDataSet, "$", tclvalue(lhsVariable), " <- factor(", .activeDataSet, "$", tclvalue(lhsVariable), ", levels=c('", paste(effLevs,sep="", collapse="', '"), "'))", sep="")
      doItAndPrint(command)
    }
    formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
    command <- paste("multinom(", formula,
                     weights, ",", " data=", .activeDataSet, subset, ", trace=FALSE)", sep="")
    #    logger(paste(modelValue, " <- ", command, sep=""))
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("summaryMultinom(", modelValue, ")", sep=""))
    activeModel(modelValue)
    if(baseValue != ""){
      command <- paste(.activeDataSet, "$", tclvalue(lhsVariable), " <- factor(", .activeDataSet, "$", tclvalue(lhsVariable), ", levels=c('", paste(effLevsO,sep="", collapse="', '"), "'))", sep="")
      doItAndPrint(command)
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="multinom", model=TRUE, reset="resetMNLNMBU", apply="multinomialLogitModelNMBU")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", row=1, column=1, columnspan=2)
  modelFormula()
  subsetBox(model=TRUE)
  tkgrid(getFrame(xBox), sticky="w", row=2, column=1, columnspan=2)
  tkgrid(outerOperatorsFrame, sticky="w", row=3, column=1, columnspan=2)
  tkgrid(formulaFrame, sticky="w", row=4, column=1, columnspan=2)
  tkgrid(subsetFrame, sticky="w", row=5, column=1, columnspan=2)
  tkgrid(labelRcmdr(baseFrame, text=gettextRcmdr("Base level (optional):")), baseV, sticky="w")
  tkgrid(baseFrame, sticky="w", row=6, column=1, columnspan=1)
#  tkgrid(labelRcmdr(weightsFrame, text=gettextRcmdr("   Weights (optional):")), weights, sticky="w")
  tkgrid(labelRcmdr(weightsFrame, text=gettextRcmdr("")), getFrame(weightComboBox), sticky="w")
  tkgrid(weightsFrame, sticky="w", row=6, column=2, columnspan=1)
  tkgrid(buttonsFrame, sticky="w", row=7, column=1, columnspan=2)
  dialogSuffix(rows=7, columns=2, focus=lhsEntry, preventDoubleClick=TRUE)
}

resetMNLNMBU <- function(){
	putRcmdr("reset.model", TRUE)
	multinomialLogitModelNMBU()
}

################################
# Customized ordinalRegression
ordinalRegressionModelNMBU <- function(){
	defaults <- list(initial.type="logistic")
	dialog.values <- getDialog("ordinalRegressionModel", defaults)
  Library("nnet")
  initializeDialog(title=gettextRcmdr("Ordinal Regression Model"))
  .activeModel <- ActiveModel()
  .activeDataSet <- ActiveDataSet()
  currentModel <- if (!is.null(.activeModel))
    class(get(.activeModel, envir=.GlobalEnv))[1] == "polr"
  else FALSE
  if (currentModel) {
    currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
    if (currentFields$data != .activeDataSet) currentModel <- FALSE
  }
	if (isTRUE(getRcmdr("reset.model"))) {
		currentModel <- FALSE
		putRcmdr("reset.model", FALSE)
	}
  UpdateModelNumber()
  modelName <- tclVar(paste("OrdRegModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  radioButtons(name="modelType",
               buttons=c("logistic", "probit", "cloglog"),
               labels=gettextRcmdr(c("Proportional-odds logit", "Ordered probit", "Complementary log-log")),
               title=gettextRcmdr("Type of Model"))
  weightsFrame <- tkframe(top)
  weightComboBox <- variableComboBox(weightsFrame, variableList=Numeric(), 
                                     initialSelection=gettextRcmdr("<no variable selected>"),
                                     title=gettextRcmdr("Weights"))
#  weightsName <- tclVar(gettextRcmdr("<none>"))
#  weights <- ttkentry(weightsFrame, width="20", textvariable=weightsName)
  onOK <- function(){
    modelValue <- trim.blanks(tclvalue(modelName))
    weight.var <- getSelection(weightComboBox)
    weights <- if (weight.var == gettextRcmdr("<no variable selected>")) ""
    else paste(" weights=", weight.var, ", ", sep="")
#    weightsValue <- tclvalue(weightsName)
    closeDialog()
    if (!is.valid.name(modelValue)){
      errorCondition(recall=ordinalRegressionModelNMBU, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    check.empty <- gsub(" ", "", tclvalue(lhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=ordinalRegressionModelNMBU, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
      return()
    }
    check.empty <- gsub(" ", "", tclvalue(rhsVariable))
    if ("" == check.empty) {
      errorCondition(recall=ordinalRegressionModelNMBU, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
      return()
    }
    if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=get(.activeDataSet, envir=.GlobalEnv)))){
      #        if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=eval(parse(text=.activeDataSet), envir=.GlobalEnv)))){
      errorCondition(recall=ordinalRegressionModelNMBU, message=gettextRcmdr("Response variable must be a factor"))
      return()
    }
    if (is.element(modelValue, listProportionalOddsModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        ordinalRegressionModelNMBU()
        return()
      }
    }
    formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
    command <- paste("polr(", formula, ', method="', tclvalue(modelTypeVariable),
                     '",', weights, ' data=', .activeDataSet, subset, ", Hess=TRUE)", sep="")
    #    logger(paste(modelValue, " <- ", command, sep=""))
    #    assign(modelValue, justDoIt(command), envir=.GlobalEnv)
    doItAndPrint(paste(modelValue, " <- ", command, sep=""))
    doItAndPrint(paste("summaryOrdinal(", modelValue, ")", sep=""))
    activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="polr", model=TRUE, reset = "resetPOLRNMBU", apply = "ordinalRegressionModelNMBU")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", row=1, column=1, columnspan=2)
  modelFormula()
  subsetBox(model=TRUE)
  tkgrid(getFrame(xBox), sticky="w", row=2, column=1, columnspan=2)
  tkgrid(outerOperatorsFrame, sticky="w", row=3, column=1, columnspan=2)
  tkgrid(formulaFrame, sticky="w", row=4, column=1, columnspan=2)
  tkgrid(subsetFrame, sticky="w", row=5, column=1, columnspan=2)
  tkgrid(modelTypeFrame, sticky="w", row=6, column=1, columnspan=1)
#  tkgrid(labelRcmdr(weightsFrame, text=gettextRcmdr("   Weights (optional):")), weights, sticky="w")
  tkgrid(labelRcmdr(weightsFrame, text=gettextRcmdr("")), getFrame(weightComboBox), sticky="w")
  tkgrid(weightsFrame, sticky="w", row=6, column=2, columnspan=1)
  tkgrid(buttonsFrame, sticky="w", row=7, column=1, columnspan=2)
  dialogSuffix(rows=7, columns=2, focus=lhsEntry, preventDoubleClick=TRUE)
}

resetPOLRNMBU <- function(){
	putRcmdr("reset.model", TRUE)
	putDialog("ordinalRegressionModelNMBU", NULL)
	ordinalRegressionModel()
}

################################
# Customized numerical summaries
numericalSummariesNMBU <- function(){
  defaults <- list(initial.x=NULL, initial.mean="1", initial.median="0", initial.sum="0", initial.sumSq="0", initial.sd="1", initial.sdErr="0", initial.var="0", initial.cv="0",
                   initial.quantiles.variable="1", 
                   initial.quantiles="0, .25, .5, .75, 1", 
                   initial.skewness="0", initial.kurtosis="0", initial.type="2",
                   initial.group=NULL)
  dialog.values <- getDialog("numericalSummariesNMBU", defaults)
  initial.group <- dialog.values$initial.group
  initializeDialog(title=gettextRcmdr("Numerical Summaries"))
  xBox <- variableListBox(top, Numeric(), selectmode="multiple", title=gettextRcmdr("Variables (pick one or more)"),
                          initialSelection=varPosn(dialog.values$initial.x, "numeric"))
  selectFrame <- tkframe(top)
  checkBoxes(frame="checkBoxFrame", boxes=c("mean", "median", "sum", "sumSq", "sd", "sdErr", "var", "cv"), 
             initialValues=c(dialog.values$initial.mean, dialog.values$initial.median, dialog.values$initial.sum, dialog.values$initial.sumSq, dialog.values$initial.sd, dialog.values$initial.sdErr, dialog.values$initial.var, dialog.values$initial.cv), #c(1,0,0,0,0,0,0,0), 
             labels=gettextRcmdr(c("Mean", "Median", "Sum", "Sum of suqares", "Standard Deviation", "Standard Error of the Mean", "Variance", "Coefficient of Variation")))
  checkBoxes(window=selectFrame, frame="skCheckBoxFrame", boxes=c("skewness", "kurtosis"), 
             initialValues=c(dialog.values$initial.skewness, dialog.values$initial.kurtosis), 
             labels=gettextRcmdr(c("Skewness", "Kurtosis")))
  radioButtons(window=selectFrame, name="typeButtons", buttons=c("b1", "b2", "b3"), values=c("1", "2", "3"), 
               initialValue=dialog.values$initial.type,
               labels=gettextRcmdr(c("Type 1", "Type 2", "Type 3")))
  quantilesVariable <- tclVar(dialog.values$initial.quantiles.variable)
  quantilesFrame <- tkframe(top)
  quantilesCheckBox <- tkcheckbutton(quantilesFrame, variable=quantilesVariable)
  quantiles <- tclVar(dialog.values$initial.quantiles)
  quantilesEntry <- ttkentry(quantilesFrame, width="20", textvariable=quantiles)
  groupsBox(recall=numericalSummariesNMBU, label=gettextRcmdr("Summarize by:"), 
            initialLabel=gettextRcmdr("Summarize by groups"), 
            initialGroup=initial.group)
  onOK <- function(){
    x <- getSelection(xBox)
    quants <- tclvalue(quantiles)
    meanVar <- tclvalue(meanVariable)
    medianVar <- tclvalue(medianVariable)
    sumVar <- tclvalue(sumVariable)
    sumSqVar <- tclvalue(sumSqVariable)
    sdVar <- tclvalue(sdVariable)
    sdErrVar <- tclvalue(sdErrVariable)
    varVar <- tclvalue(varVariable)
    cvVar <- tclvalue(cvVariable)
    quantsVar <- tclvalue(quantilesVariable)
    skewnessVar <- tclvalue(skewnessVariable)
    kurtosisVar <- tclvalue(kurtosisVariable)
    typeVar <- tclvalue(typeButtonsVariable)
    putDialog("numericalSummariesNMBU", list(
      initial.x=x, initial.mean=meanVar, initial.median=medianVar, initial.sum=sumVar, initial.sumSq=sumSqVar, initial.sd=sdVar, initial.sdErr=sdErrVar, initial.var=sdVar, initial.cv=cvVar,
      initial.quantiles.variable=quantsVar, initial.quantiles=quants,
      initial.skewness=skewnessVar, initial.kurtosis=kurtosisVar, initial.type=typeVar,
      initial.group=if (.groups != FALSE) .groups else NULL
    ))		
    if (length(x) == 0){
      errorCondition(recall=numericalSummariesNMBU, message=gettextRcmdr("You must select a variable."))
      return()
    }
    closeDialog()
    quants <- paste("c(", gsub(",+", ",", gsub(" ", ",", quants)), ")", sep="")
    .activeDataSet <- ActiveDataSet()
    vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
    else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
    vars <- paste(.activeDataSet, "[,", vars, "]", sep="")
    stats <- paste("c(",
                   paste(c('"mean"', '"median"', '"sum"', '"sumSq"', '"sd"', '"sdErr"', '"var"', '"quantiles"', '"cv"', '"skewness"', '"kurtosis"')
                         [c(meanVar, medianVar, sumVar, sumSqVar, sdVar, sdErrVar, varVar, quantsVar, cvVar, skewnessVar, kurtosisVar) == 1], 
                         collapse=", "), ")", sep="")
    if (stats == "c()"){
      errorCondition(recall=numericalSummariesNMBU, message=gettextRcmdr("No statistics selected."))
      return()
    }
    type.text <- if (skewnessVar == 1 || kurtosisVar == 1 || quantsVar == 1) paste(', type="', typeVar, '"', sep="") else ""
    command <- if (.groups != FALSE) {
      grps <- paste(.activeDataSet, "$", .groups, sep="")
      paste("numSummaryNMBU(", vars, ", groups=", grps, ", statistics=", stats, 
            ", quantiles=", quants, type.text, ")", sep="")
    }
    else  paste("numSummaryNMBU(", vars, ", statistics=", stats, 
                ", quantiles=", quants, type.text, ")", sep="")
    doItAndPrint(command) 
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="numSummary", reset="numericalSummariesNMBU")
  tkgrid(getFrame(xBox), sticky="nw")    
  tkgrid(checkBoxFrame, sticky="w")
  tkgrid(skCheckBoxFrame, typeButtonsFrame, sticky="nw")
  tkgrid(selectFrame, sticky="w")
  tkgrid(labelRcmdr(quantilesFrame, text=gettextRcmdr("Quantiles")), quantilesCheckBox,
         labelRcmdr(quantilesFrame, text=gettextRcmdr(" quantiles:")), quantilesEntry, sticky="w")
  tkgrid(quantilesFrame, sticky="w")
  tkgrid(groupsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=7, columns=1)
}

#####################################
# Bugfixed one-way ANOVA
oneWayAnovaNMBU <- function(){
  initializeDialog(title=gettextRcmdr("One-Way Analysis of Variance"))
  UpdateModelNumber()
  modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  groupBox <- variableListBox(top, Factors(), title=gettextRcmdr("Groups (pick one)"))
  responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
  optionsFrame <- tkframe(top)
  pairwiseVariable <- tclVar("0")
  pairwiseCheckBox <- tkcheckbutton(optionsFrame, variable=pairwiseVariable)
  onOK <- function(){
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      UpdateModelNumber(-1)
      errorCondition(recall=oneWayAnova, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listAOVModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        tkdestroy(top)
        oneWayAnova()
        return()
      }
    }
    putRcmdr("modelWithSubset", FALSE)
    group <- getSelection(groupBox)
    response <- getSelection(responseBox)
    closeDialog()
    if (length(group) == 0){
      errorCondition(recall=oneWayAnovaNMBU, message=gettextRcmdr("You must select a groups factor."))
      return()
    }
    if (length(response) == 0){
      errorCondition(recall=oneWayAnovaNMBU, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    command <- paste(modelValue, " <- aov(", response, " ~ ", group, ", data=", .activeDataSet, ")", sep="")
    justDoIt(command)
    logger(command)
    doItAndPrint(paste("summary(", modelValue, ")", sep=""))
    doItAndPrint(paste("numSummary(", .activeDataSet, "$", response, " , groups=", .activeDataSet, "$", group,
                       ', statistics=c("mean", "sd"))', sep=""))
    activeModel(modelValue)
    pairwise <- tclvalue(pairwiseVariable)
    if (pairwise == 1) {
      if (eval(parse(text=paste("length(levels(", .activeDataSet, "$", group, ")) < 3"))))
        Message(message=gettextRcmdr("Factor has fewer than 3 levels; pairwise comparisons omitted."),
                type="warning")
      # the following lines modified by Richard Heiberger and subsequently by J. Fox
      else {
        command <- paste(".Pairs <- glht(", modelValue, ", linfct = mcp(", group, ' = "Tukey"))', sep="")
        justDoIt(command)
        logger(command)
        doItAndPrint("summary(.Pairs) # pairwise tests")
        doItAndPrint("confint(.Pairs) # confidence intervals")
        doItAndPrint("cld(.Pairs) # compact letter display")
        justDoIt("old.oma <- par(oma=c(0,5,0,0))")
        logger("old.oma <- par(oma=c(0,5,0,0))")
        justDoIt("plot(confint(.Pairs))")
        logger("plot(confint(.Pairs))")
        justDoIt("par(old.oma)")
        logger("par(old.oma)")
        logger("remove(.Pairs)")
        remove(.Pairs, envir=.GlobalEnv)
      }
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="anova", model=TRUE)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", columnspan=2)
  tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
  tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Pairwise comparisons of means")), pairwiseCheckBox, sticky="w")
  tkgrid(optionsFrame, sticky="w", columnspan=2)
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}

#####################################
# Bugfixed multi-way ANOVA
multiWayAnovaNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Multi-Way Analysis of Variance"))
  UpdateModelNumber()
  modelName <- tclVar(paste("AnovaModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  groupBox <- variableListBox(top, Factors(), selectmode="multiple", title=gettextRcmdr("Factors (pick one or more)"))
  responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
  onOK <- function(){
    modelValue <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(modelValue)){
      UpdateModelNumber(-1)
      errorCondition(recall=multiWayAnovaNMBU, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
      return()
    }
    if (is.element(modelValue, listAOVModels())) {
      if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        tkdestroy(top)
        multiWayAnova()
        return()
      }
    }
    putRcmdr("modelWithSubset", FALSE)
    groups <- getSelection(groupBox)
    response <- getSelection(responseBox)
    closeDialog()
    if (length(groups) == 0){
      errorCondition(recall=multiWayAnovaNMBU, message=gettextRcmdr("You must select at least one factor."))
      return()
    }
    if (length(response) == 0){
      errorCondition(recall=multiWayAnovaNMBU, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    groups.list <- paste(paste(groups, "=", .activeDataSet, "$", groups, sep=""), collapse=", ")
    doItAndPrint(paste(modelValue, " <- (lm(", response, " ~ ", paste(groups, collapse="*"),
                       ", data=", .activeDataSet, "))", sep=""))
    doItAndPrint(paste("Anova(", modelValue, ")", sep=""))
    doItAndPrint(paste("tapply(", .activeDataSet, "$", response, ", list(", groups.list,
                       "), mean, na.rm=TRUE) # means", sep=""))
    doItAndPrint(paste("tapply(", .activeDataSet, "$", response, ", list(", groups.list,
                       "), sd, na.rm=TRUE) # std. deviations", sep=""))
    doItAndPrint(paste("tapply(", .activeDataSet, "$", response, ", list(", groups.list,
                       "), function(x) sum(!is.na(x))) # counts", sep=""))
    activeModel(modelValue)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="Anova", model=TRUE)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model: ")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", columnspan=2)
  tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}

#################################
# Two sample t-test (summarized data)
twoSamplesTTestSum <- function(){
  initializeDialog(title=gettextRcmdr("Two samples t-test"))
  onOK <- function(){ # Actions to perform
    mu1 <- as.character(tclvalue(mu1Level))
    n1  <- as.character(tclvalue(n1Level))
    sd1 <- as.character(tclvalue(sd1Level))
    mu2 <- as.character(tclvalue(mu2Level))
    n2  <- as.character(tclvalue(n2Level))
    sd2 <- as.character(tclvalue(sd2Level))
    mu0 <- as.character(tclvalue(mu0Level))
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    closeDialog()
    doItAndPrint(paste("t_test_sum(means=c(", mu1,",",mu2, "), sds=c(", sd1,",",sd2, "), ns=c(", n1,",",n2, "), alternative='",
                       alternative, "', mu=", mu0, ", var.equal=", variances, ", conf.level=", level, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="t_test")
  
  grFrame <- tkframe(top);
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr(""), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Group 1"), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Group 2"), fg="blue"),sticky="w")
  tkgrid(grFrame, sticky="nw", row=1, column=1)
  
  muFrame <- tkframe(top)
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Mean"), fg="blue"),sticky="w")
  mu1Level <- tclVar("");   mu1Field <- ttkentry(muFrame, width="6", textvariable=mu1Level)
  mu2Level <- tclVar("");   mu2Field <- ttkentry(muFrame, width="6", textvariable=mu2Level)
  tkgrid(mu1Field, sticky="nw");  tkgrid(mu2Field, sticky="nw");  tkgrid(muFrame, sticky="nw", row=1, column=2)
  
  nFrame <- tkframe(top)
  tkgrid(labelRcmdr(nFrame, text=gettextRcmdr("N"), fg="blue"),sticky="w")
  n1Level <- tclVar("");   n1Field <- ttkentry(nFrame, width="6", textvariable=n1Level)
  n2Level <- tclVar("");   n2Field <- ttkentry(nFrame, width="6", textvariable=n2Level)
  tkgrid(n1Field, sticky="nw");  tkgrid(n2Field, sticky="nw");  tkgrid(nFrame, sticky="nw", row=1, column=3)
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("St.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  sd2Level <- tclVar("");   sd2Field <- ttkentry(sdFrame, width="6", textvariable=sd2Level)
  tkgrid(sd1Field, sticky="nw");  tkgrid(sd2Field, sticky="nw");  tkgrid(sdFrame, sticky="nw", row=1, column=4)
  
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu1-mu2"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw", row=2, column=1)
  
  
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu1-mu2", "Difference > mu1-mu2")), title=gettextRcmdr("Alternative hypothesis"))
  tkgrid(alternativeFrame, sticky="nw",row=3, column=1, columnspan=3)
  
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence level   "), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="nw",row=3, column=4, columnspan=2)
  radioButtons(top, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(variancesFrame, sticky="nw",row=3, column=6, columnspan=2)
  tkgrid(buttonsFrame, sticky="w", row=4,column=1, columnspan=6)
  dialogSuffix(rows=4, columns=7)
}


#################################
# One sample t-test (summarized data)
oneSamplesTTestSum <- function(){
  initializeDialog(title=gettextRcmdr("One sample t-test"))
  onOK <- function(){ # Actions to perform
    mu1 <- as.character(tclvalue(mu1Level))
    n1  <- as.character(tclvalue(n1Level))
    sd1 <- as.character(tclvalue(sd1Level))
    mu0 <- as.character(tclvalue(mu0Level))
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    closeDialog()
    doItAndPrint(paste("t_test_sum(means=", mu1, ", sds=", sd1, ", ns=", n1, ", alternative='",
                       alternative, "', mu=", mu0, ", conf.level=", level, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="t_test")
  
  grFrame <- tkframe(top);
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr(""), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Summaries"), fg="blue"),sticky="w")
  tkgrid(grFrame, sticky="nw", row=1, column=1)
  
  muFrame <- tkframe(top)
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Mean"), fg="blue"),sticky="w")
  mu1Level <- tclVar("");   mu1Field <- ttkentry(muFrame, width="6", textvariable=mu1Level)
  tkgrid(mu1Field, sticky="nw");  tkgrid(muFrame, sticky="nw", row=1, column=2)
  
  nFrame <- tkframe(top)
  tkgrid(labelRcmdr(nFrame, text=gettextRcmdr("N"), fg="blue"),sticky="w")
  n1Level <- tclVar("");   n1Field <- ttkentry(nFrame, width="6", textvariable=n1Level)
  tkgrid(n1Field, sticky="nw");  tkgrid(nFrame, sticky="nw", row=1, column=3)
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("St.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  tkgrid(sd1Field, sticky="nw");  tkgrid(sdFrame, sticky="nw", row=1, column=4)
  
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu0"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw", row=2, column=1)
  
  
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "True mean < mu0", "True mean > mu0")), title=gettextRcmdr("Alternative hypothesis"))
  tkgrid(alternativeFrame, sticky="nw",row=3, column=1, columnspan=3)
  
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="nw",row=3, column=4, columnspan=2)
  tkgrid(buttonsFrame, sticky="w", row=4,column=1, columnspan=6)
  dialogSuffix(rows=4, columns=7)
}

#################################
# Independent samples t-test
independentSamplesTTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Independent Samples t-Test"))
  variablesFrame <- tkframe(top)
  groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), title=gettextRcmdr("Groups (pick one)"))
  responseBox <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
  onOK <- function(){
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall=independentSamplesTTestNMBU, message=gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall=independentSamplesTTestNMBU, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    mu0 <- as.character(tclvalue(mu0Level))
    closeDialog()
    doItAndPrint(paste("t_test(", response, "~", group,
                       ", alternative='", alternative, "', mu=", mu0, ", conf.level=", level,
                       ", var.equal=", variances,
                       ", data=", ActiveDataSet(), ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_test")
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  radioButtons(optionsFrame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(responseBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu1-mu2"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  groupsLabel(groupsBox=groupBox)
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), confidenceFrame, labelRcmdr(optionsFrame, text="    "),
         variancesFrame, sticky="nw")
  tkgrid(optionsFrame, sticky="nw")
  
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=6, columns=1)
}

#################################
# Paired samples t-test
pairedTTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Paired t-Test"))
  .numeric <- Numeric()
  xBox <- variableListBox(top, .numeric, title=gettextRcmdr("First variable (pick one)"))
  yBox <- variableListBox(top, .numeric, title=gettextRcmdr("Second variable (pick one)"))
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    if (length(x) == 0 | length(y) == 0){
      errorCondition(recall=pairedTTestNMBU, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y){
      errorCondition(recall=pairedTTestNMBU, message=gettextRcmdr("Variables must be different."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    mu <- tclvalue(muVariable)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("t_test(", .activeDataSet, "$", x, ", ",
                       .activeDataSet, "$", y,
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
                       ", paired=TRUE)", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_test")
  muFrame <- tkframe(top)
  muVariable <- tclVar("0.0")
  muField <- ttkentry(muFrame, width="8", textvariable=muVariable)
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Null hypothesis: mu = ")), muField, sticky="w")
  tkgrid(muFrame, sticky="w")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"))
  tkgrid(confidenceField, sticky="w")
  tkgrid(alternativeFrame, confidenceFrame, sticky="nw")
  
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=4, columns=2)
}

#################################
# Single sample t-test
singleSampleTTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Single-Sample t-Test"))
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  onOK <- function(){
    x <- getSelection(xBox)
    if (length(x) == 0){
      errorCondition(recall=singleSampleTTestNMBU, message=gettextRcmdr("You must select a variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    mu <- tclvalue(muVariable)
    closeDialog()
    doItAndPrint(paste("t_test(", ActiveDataSet(), "$", x,
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
                       ")", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_test")
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Population mean != mu0", "Population mean < mu0", "Population mean > mu0")),
               title=gettextRcmdr("Alternative Hypothesis"))
  rightFrame <- tkframe(top)
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  muFrame <- tkframe(rightFrame)
  muVariable <- tclVar("0.0")
  muField <- ttkentry(muFrame, width="8", textvariable=muVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelRcmdr(rightFrame, text=""), sticky="w")
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Null hypothesis: mu = ")), muField, sticky="w")
  tkgrid(muFrame, sticky="w")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level: ")), confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="w")
  tkgrid(alternativeFrame, rightFrame, sticky="nw")
  
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(confidenceField, sticky="e")
  dialogSuffix(rows=4, columns=2)
}

#################################
# Two sample z-test (summarized data)
twoSamplesZTestSum <- function(){
  initializeDialog(title=gettextRcmdr("Two samples z-test"))
  onOK <- function(){ # Actions to perform
    mu1 <- as.character(tclvalue(mu1Level))
    n1  <- as.character(tclvalue(n1Level))
    sd1 <- as.character(tclvalue(sd1Level))
    mu2 <- as.character(tclvalue(mu2Level))
    n2  <- as.character(tclvalue(n2Level))
    sd2 <- as.character(tclvalue(sd2Level))
    mu0 <- as.character(tclvalue(mu0Level))
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    closeDialog()
    doItAndPrint(paste("z_test_sum(means=c(", mu1,",",mu2, "), sds=c(", sd1,",",sd2, "), ns=c(", n1,",",n2, "), alternative='",
                       alternative, "', mu=", mu0, ", var.equal=", variances, ", conf.level=", level, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="z_test")
  
  grFrame <- tkframe(top);
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr(""), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Group 1"), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Group 2"), fg="blue"),sticky="w")
  tkgrid(grFrame, sticky="nw", row=1, column=1)
  
  muFrame <- tkframe(top)
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Mean"), fg="blue"),sticky="w")
  mu1Level <- tclVar("");   mu1Field <- ttkentry(muFrame, width="6", textvariable=mu1Level)
  mu2Level <- tclVar("");   mu2Field <- ttkentry(muFrame, width="6", textvariable=mu2Level)
  tkgrid(mu1Field, sticky="nw");  tkgrid(mu2Field, sticky="nw");  tkgrid(muFrame, sticky="nw", row=1, column=2)
  
  nFrame <- tkframe(top)
  tkgrid(labelRcmdr(nFrame, text=gettextRcmdr("N"), fg="blue"),sticky="w")
  n1Level <- tclVar("");   n1Field <- ttkentry(nFrame, width="6", textvariable=n1Level)
  n2Level <- tclVar("");   n2Field <- ttkentry(nFrame, width="6", textvariable=n2Level)
  tkgrid(n1Field, sticky="nw");  tkgrid(n2Field, sticky="nw");  tkgrid(nFrame, sticky="nw", row=1, column=3)
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("St.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  sd2Level <- tclVar("");   sd2Field <- ttkentry(sdFrame, width="6", textvariable=sd2Level)
  tkgrid(sd1Field, sticky="nw");  tkgrid(sd2Field, sticky="nw");  tkgrid(sdFrame, sticky="nw", row=1, column=4)
  
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu1-mu2"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw", row=2, column=1)
  
  
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu1-mu2", "Difference > mu1-mu2")), title=gettextRcmdr("Alternative hypothesis"))
  tkgrid(alternativeFrame, sticky="nw",row=3, column=1, columnspan=3)
  
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence level   "), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="nw",row=3, column=4, columnspan=2)
  radioButtons(top, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(variancesFrame, sticky="nw",row=3, column=6, columnspan=2)
  tkgrid(buttonsFrame, sticky="w", row=4,column=1, columnspan=6)
  dialogSuffix(rows=4, columns=7)
}


#################################
# One sample z-test (summarized data)
oneSamplesZTestSum <- function(){
  initializeDialog(title=gettextRcmdr("One sample z-test"))
  onOK <- function(){ # Actions to perform
    mu1 <- as.character(tclvalue(mu1Level))
    n1  <- as.character(tclvalue(n1Level))
    sd1 <- as.character(tclvalue(sd1Level))
    mu0 <- as.character(tclvalue(mu0Level))
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    closeDialog()
    doItAndPrint(paste("z_test_sum(means=", mu1, ", sds=", sd1, ", ns=", n1, ", alternative='",
                       alternative, "', mu=", mu0, ", conf.level=", level, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="z_test")
  
  grFrame <- tkframe(top);
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr(""), fg="blue"),sticky="w")
  tkgrid(labelRcmdr(grFrame, text=gettextRcmdr("Summaries"), fg="blue"),sticky="w")
  tkgrid(grFrame, sticky="nw", row=1, column=1)
  
  muFrame <- tkframe(top)
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Mean"), fg="blue"),sticky="w")
  mu1Level <- tclVar("");   mu1Field <- ttkentry(muFrame, width="6", textvariable=mu1Level)
  tkgrid(mu1Field, sticky="nw");  tkgrid(muFrame, sticky="nw", row=1, column=2)
  
  nFrame <- tkframe(top)
  tkgrid(labelRcmdr(nFrame, text=gettextRcmdr("N"), fg="blue"),sticky="w")
  n1Level <- tclVar("");   n1Field <- ttkentry(nFrame, width="6", textvariable=n1Level)
  tkgrid(n1Field, sticky="nw");  tkgrid(nFrame, sticky="nw", row=1, column=3)
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("St.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  tkgrid(sd1Field, sticky="nw");  tkgrid(sdFrame, sticky="nw", row=1, column=4)
  
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu0"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw", row=2, column=1)
  
  
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "True mean < mu0", "True mean > mu0")), title=gettextRcmdr("Alternative hypothesis"))
  tkgrid(alternativeFrame, sticky="nw",row=3, column=1, columnspan=3)
  
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="nw",row=3, column=4, columnspan=2)
  tkgrid(buttonsFrame, sticky="w", row=4,column=1, columnspan=6)
  dialogSuffix(rows=4, columns=7)
}

#################################
# Independent samples z-test
independentSamplesZTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Independent Samples z-Test"))
  variablesFrame <- tkframe(top)
  groupBox <- variableListBox(variablesFrame, TwoLevelFactors(), title=gettextRcmdr("Groups (pick one)"))
  responseBox <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
  onOK <- function(){
    sd1 <- as.character(tclvalue(sd1Level))
    sd2 <- as.character(tclvalue(sd2Level))
    if(sd1==gettextRcmdr("") || sd1==gettextRcmdr("")){
      errorCondition(recall=independentSamplesZTestNMBU, message=gettextRcmdr("You must specify both standard deviations."))
    }
    sd1 <- as.numeric(sd1)
    sd2 <- as.numeric(sd2)
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall=independentSamplesZTestNMBU, message=gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall=independentSamplesZTestNMBU, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    mu0 <- as.character(tclvalue(mu0Level))
    closeDialog()
    doItAndPrint(paste("z_test(", response, "~", group,
                       ", alternative='", alternative, "', mu=", mu0, ", conf.level=", level,
                       ", var.equal=", variances,
                       ", data=", ActiveDataSet(), ", sds=c(", sd1, ",", sd2, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_test")
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  radioButtons(optionsFrame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(responseBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu1-mu2"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  groupsLabel(groupsBox=groupBox)
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), confidenceFrame, labelRcmdr(optionsFrame, text="    "),
         variancesFrame, sticky="nw")
  tkgrid(optionsFrame, sticky="nw")
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("Known st.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  sd2Level <- tclVar("");   sd2Field <- ttkentry(sdFrame, width="6", textvariable=sd2Level)
  tkgrid(sd1Field, sd2Field, sticky="w");  tkgrid(sdFrame, sticky="nw")
  
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=6, columns=1)
}

#################################
# Paired samples z-test
pairedZTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Paired z-Test"))
  .numeric <- Numeric()
  xBox <- variableListBox(top, .numeric, title=gettextRcmdr("First variable (pick one)"))
  yBox <- variableListBox(top, .numeric, title=gettextRcmdr("Second variable (pick one)"))
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    sd1 <- as.character(tclvalue(sd1Level))
    z.test <- FALSE
    if(sd1==gettextRcmdr("")){
      errorCondition(recall=pairedZTestNMBU, message=gettextRcmdr("You must specify a standard deviation."))
    }
    sd1 <- as.numeric(sd1)
    if (length(x) == 0 | length(y) == 0){
      errorCondition(recall=pairedZTestNMBU, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y){
      errorCondition(recall=pairedZTestNMBU, message=gettextRcmdr("Variables must be different."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    mu <- tclvalue(muVariable)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("z_test(", .activeDataSet, "$", x, ", ",
                       .activeDataSet, "$", y,
                       ", alternative='", alternative, "', conf.level=", level,
                       ", paired=TRUE, sds=", sd1, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_test")
  muFrame <- tkframe(top)
  muVariable <- tclVar("0.0")
  muField <- ttkentry(muFrame, width="8", textvariable=muVariable)
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(top)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Null hypothesis: mu = ")), muField, sticky="w")
  tkgrid(muFrame, sticky="w")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"))
  tkgrid(confidenceField, sticky="w")
  tkgrid(alternativeFrame, confidenceFrame, sticky="nw")
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("Known st.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  tkgrid(sd1Field, sticky="w");  tkgrid(sdFrame, sticky="nw")
  
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=5, columns=2)
}

#################################
# Single sample z-test
singleSampleZTestNMBU <- function(){
  initializeDialog(title=gettextRcmdr("Single-Sample z-Test"))
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
  onOK <- function(){
    x <- getSelection(xBox)
    sd1 <- as.character(tclvalue(sd1Level))
    if(sd1==gettextRcmdr("")){
      errorCondition(recall=singleSampleZTestNMBU, message=gettextRcmdr("You must specify a standard deviation."))
    }
    sd1 <- as.numeric(sd1)
    if (length(x) == 0){
      errorCondition(recall=singleSampleZTestNMBU, message=gettextRcmdr("You must select a variable."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    mu <- tclvalue(muVariable)
    closeDialog()
    doItAndPrint(paste("z_test(", ActiveDataSet(), "$", x,
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
                       ", sds=", sd1,")", sep=""))
    tkdestroy(top)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="z_test")
  radioButtons(top, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Population mean != mu0", "Population mean < mu0", "Population mean > mu0")),
               title=gettextRcmdr("Alternative Hypothesis"))
  rightFrame <- tkframe(top)
  confidenceFrame <- tkframe(rightFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  muFrame <- tkframe(rightFrame)
  muVariable <- tclVar("0.0")
  muField <- ttkentry(muFrame, width="8", textvariable=muVariable)
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelRcmdr(rightFrame, text=""), sticky="w")
  tkgrid(labelRcmdr(muFrame, text=gettextRcmdr("Null hypothesis: mu = ")), muField, sticky="w")
  tkgrid(muFrame, sticky="w")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level: ")), confidenceField, sticky="w")
  tkgrid(confidenceFrame, sticky="w")
  tkgrid(alternativeFrame, rightFrame, sticky="nw")
  
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("Known st.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  tkgrid(sd1Field, sticky="w");  tkgrid(sdFrame, sticky="nw")
  
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid.configure(confidenceField, sticky="e")
  dialogSuffix(rows=5, columns=2)
}


#################################
# Unstacked two sample t-test
twoSamplesTTest <- function(){
  initializeDialog(title=gettextRcmdr("Two Samples t-Test"))
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("First variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Second variable (pick one)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    if (length(x) == 0 | length(y) == 0){
      errorCondition(recall=twoSamplesTTest, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y){
      errorCondition(recall=twoSamplesTTest, message=gettextRcmdr("Variables must be different."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    mu <- tclvalue(mu0Level)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("t_test(", .activeDataSet, "$", x, ", ",
                       .activeDataSet, "$", y,
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
                       ", var.equal=", variances,
                       ", data=", ActiveDataSet(), ")", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="t_test")
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  radioButtons(optionsFrame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu0"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), confidenceFrame, labelRcmdr(optionsFrame, text="    "),
         variancesFrame, sticky="nw")
  tkgrid(optionsFrame, sticky="nw")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=4, columns=1)
}

#################################
# Unstacked two sample z-test
twoSamplesZTest <- function(){
  initializeDialog(title=gettextRcmdr("Two Samples z-Test"))
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("First variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Second variable (pick one)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    sd1 <- as.character(tclvalue(sd1Level))
    sd2 <- as.character(tclvalue(sd2Level))
    if(sd1==gettextRcmdr("") || sd1==gettextRcmdr("")){
      errorCondition(recall=twoSamplesZTest, message=gettextRcmdr("You must specify both standard deviations."))
    }
    sd1 <- as.numeric(sd1)
    sd2 <- as.numeric(sd2)
    if (length(x) == 0 | length(y) == 0){
      errorCondition(recall=twoSamplesZTest, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y){
      errorCondition(recall=twoSamplesZTest, message=gettextRcmdr("Variables must be different."))
      return()
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    level <- tclvalue(confidenceLevel)
    variances <- as.character(tclvalue(variancesVariable))
    mu <- tclvalue(mu0Level)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("z_test(", .activeDataSet, "$", x, ", ",
                       .activeDataSet, "$", y,
                       ", alternative='", alternative, "', mu=", mu, ", conf.level=", level,
                       ", var.equal=", variances,
                       ", data=", ActiveDataSet(), ", sds=c(", sd1, ",", sd2, "))", sep=""))
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="z_test")
  optionsFrame <- tkframe(top)
  radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < mu0", "Difference > mu0")), title=gettextRcmdr("Alternative Hypothesis"))
  confidenceFrame <- tkframe(optionsFrame)
  confidenceLevel <- tclVar(".95")
  confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel)
  radioButtons(optionsFrame, name="variances", buttons=c("yes", "no"), values=c("TRUE", "FALSE"), initialValue="FALSE",
               labels=gettextRcmdr(c("Yes", "No")), title=gettextRcmdr("Assume equal variances?"))
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  mu0Frame <- tkframe(top)
  tkgrid(labelRcmdr(mu0Frame, text=gettextRcmdr("mu0"), fg="blue"),sticky="w")
  mu0Level <- tclVar("0");    mu0Field <- ttkentry(mu0Frame, width="6", textvariable=mu0Level)
  tkgrid(mu0Field, sticky="nw");  tkgrid(mu0Frame, sticky="nw")
  tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w")
  tkgrid(confidenceField, sticky="w")
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), confidenceFrame, labelRcmdr(optionsFrame, text="    "),
         variancesFrame, sticky="nw")
  tkgrid(optionsFrame, sticky="nw")
  sdFrame <- tkframe(top)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("Known st.dev."), fg="blue"),sticky="w")
  sd1Level <- tclVar("");   sd1Field <- ttkentry(sdFrame, width="6", textvariable=sd1Level)
  sd2Level <- tclVar("");   sd2Field <- ttkentry(sdFrame, width="6", textvariable=sd2Level)
  tkgrid(sd1Field, sd2Field, sticky="w");  tkgrid(sdFrame, sticky="nw")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=5, columns=1)
}


#####################################
# Create simplex mixture design
simplex.analysis <- function(){
  initializeDialog(title=gettextRcmdr("Create simplex mixture design"))
  .numeric <- Numeric()
  plotFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  x1Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Left variable (pick one)"))
  x2Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Top variable (pick one)"))
  x3Box <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Right variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Response variable (pick one)"))
  UpdateModelNumber()
  modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  formatFrame <- tkframe(top)
  
  onOK <- function(){ # Actions to perform
    x1 <- getSelection(x1Box)
    x2 <- getSelection(x2Box)
    x3 <- getSelection(x3Box)    
    y <- getSelection(yBox)
    do.plot <- tclvalue(plotVariable)
    linear <- as.character(tclvalue(linearVariable))
    modelValue <- trim.blanks(tclvalue(modelName))
    closeDialog()
    if (0 == length(x1) || 0 == length(x2) || 0 == length(x3) || 0 == length(y)) {
      errorCondition(recall=simplex.analysis, message=gettextRcmdr("Variables must be chosen from all boxes."))
      return()
    }
    if (!is.valid.name(modelValue)){
      errorCondition(recall=simplex.analysis, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      return()
    }
    .activeDataSet <- ActiveDataSet()
    if(linear==gettextRcmdr("linear")){
      formula1 <- paste("lm(", paste("formula(",y," ~ ",x1," + ",x2," + ",x3, " -1)", sep=""), sep="")
      command <- paste(formula1, ", data=", .activeDataSet,")", sep="")
      #      logger(paste(modelValue, " <- ", command, sep=""))
      #      assign(modelValue, justDoIt(command), envir=.GlobalEnv)
      doItAndPrint(paste(modelValue, " <- ", command, sep=""))
      activeModel(modelValue)
      doItAndPrint(paste("Anova(", modelValue, ", type='II')", sep=""))
    } else {
      formula1 <- paste("lm(", paste("formula(",y," ~ (",x1," + ",x2," + ",x3, ")^2 -1)", sep=""), sep="")
      command <- paste(formula1, ", data=", .activeDataSet,")", sep="")
      #      logger(paste(modelValue, " <- ", command, sep=""))
      #      assign(modelValue, justDoIt(command), envir=.GlobalEnv)
      doItAndPrint(paste(modelValue, " <- ", command, sep=""))
      activeModel(modelValue)
      doItAndPrint(paste("Anova(", modelValue, ", type='II')", sep=""))
    }
    if(do.plot == gettextRcmdr("1")){
      mixtureGUI()}		
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="plot")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w", column=1, row=1, columnspan=2)
  tkgrid(getFrame(x1Box), labelRcmdr(variablesFrame, text="    "), getFrame(x2Box), labelRcmdr(variablesFrame, text="    "), getFrame(x3Box), sticky="nw")
  tkgrid(getFrame(yBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw", row=2, column=1, columnspan=2)
  radioButtonsNMBU(formatFrame,name="linear", buttons=c("linear", "quadratic"), values=c("linear", "quadratic"), initialValue = "linear",
                   labels=gettextRcmdr(c("Linear model", "Quadratic model")))
  tkgrid(formatFrame, row=3, column=1, columnspan=1, rowspan=1, sticky="w")
  checkBoxes(frame="plotFrame", boxes=c("plot"), initialValues=c("0"), labels=gettextRcmdr(c("Plot responce surface")))
  tkgrid(plotFrame, row=4, column=1, columnspan=1, rowspan=1, sticky="w")
  tkgrid(buttonsFrame, sticky="w", row=5, column=1, columnspan=1)
  dialogSuffix(rows=5, columns=1)
}



#####################################
# Relevant Component Analysis (after Helland and Almøy, 1994)
RelComp <- function(){
  initializeDialog(title=gettextRcmdr("Relevant Components Plot"))
  variablesFrame1 <- tkframe(top)
  .numeric <- Numeric()
  .variable <- Variables()
  xBox <- variableListBox(variablesFrame1, .numeric, selectmode="multiple",title=gettextRcmdr("Explanatory variables (pick two or more)"))
  yBox <- variableListBox(variablesFrame1, .variable, title=gettextRcmdr("Response variable (pick one)"))
  subsetBox()
  compFrame <- tkframe(top)  
  compVar <- tclVar("3")
  compEntry <- ttkentry(compFrame, width="3", textvariable=compVar)
  
  checkBoxes(frame="optionsFrame", boxes=c("center", "scale"), initialValues=c("1", "0"),
             labels=gettextRcmdr(c("Center predictors", "Scale predictors")))
  
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    nvar <- length(x)
    y <- getSelection(yBox)
    center <- tclvalue(centerVariable)
    scale <- tclvalue(scaleVariable)
    subset <- tclvalue(subsetVariable)
    ncomp <- tclvalue(compVar)
    closeDialog()
    if (2 > nvar) {
      errorCondition(recall=Relcomp, message=gettextRcmdr("Fewer than 2 variables selected."))
      return()
    }    
    if (is.element(y, x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=Relcomp, message=gettextRcmdr("Response and explanatory variables must be different."))
      return()
    }
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ", subset=NULL" else paste(", subset=", subset, sep="")
    if(trim.blanks(ncomp) == gettextRcmdr("")){
      ncomp <- nvar
    }
    ncomp <- paste(", ncomp=",ncomp,sep="")
    
    .activeDataSet <- ActiveDataSet()
    docenter <- ifelse(center == "1",TRUE, FALSE)
    doscale <- ifelse(scale == "1",TRUE, FALSE)
    
    x <- paste('"', x, '"', sep="")  
    getY <- paste(.activeDataSet,"$",y,sep="")
    getX <- paste(.activeDataSet, "[,c(",paste(x, collapse=","),")]",sep="")
    command <- paste("plotprops(",getY,",",getX,", doscaleX=",doscale,",docenterX=",docenter,ncomp,subset,")",sep="")  
    justDoIt(command)
    logger(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="plot", model=TRUE)
  tkgrid(variablesFrame1, row=1, columnspan=2, sticky="n")
  tkgrid(getFrame(yBox), row=1, column=1, sticky="nw")
  tkgrid(getFrame(xBox), row=1, column=2, sticky="nw")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(labelRcmdr(compFrame, text=gettextRcmdr("Number of components")), compEntry, sticky="w")
  tkgrid(compFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=6, columns=2)  
}


#####################################
# Tally of discrete variable
tally.GUI <- function(){
  initializeDialog(title=gettextRcmdr("Tally of discrete variable"))
  .numeric <- Numeric()
  plotFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Tally variable (pick one)"))
  
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    if (1 != length(x)) {
      errorCondition(recall=tally.GUI, message=gettextRcmdr("Exactly one variable must be chosen"))
      return()
    }
    command <- paste("tally(",ActiveDataSet(), "$", x,")", sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="summary")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw", row=1, column=1, columnspan=1)
  tkgrid(buttonsFrame, sticky="w", row=2, column=1, columnspan=1)
  dialogSuffix(rows=2, columns=1)
}

####################################
# Power computation for t statistics
powerTtest <- function(){
  initializeDialog(title=gettextRcmdr("Power calculations for one and two sample t tests"))
  onOK <- function(){ # Actions to perform
    nVal <- tclvalue(nLevel)
    deltaVal <- tclvalue(deltaLevel)
    sdVal <- tclvalue(sdLevel)
    sigVal <- tclvalue(sigLevel)
    powVal <- tclvalue(powLevel)
    nBlank <- 0
    if(trim.blanks(nVal) != gettextRcmdr("")){
      nBlank <- nBlank + 1
    } else {nVal <- "NULL"}
    if(trim.blanks(deltaVal) != gettextRcmdr("")){
      nBlank <- nBlank + 1
    } else {deltaVal <- "NULL"}
    if(trim.blanks(sdVal) != gettextRcmdr("")){
      nBlank <- nBlank + 1
    } else {sdVal <- "NULL"}
    if(trim.blanks(sigVal) != gettextRcmdr("")){
      nBlank <- nBlank + 1
    } else {sigVal <- "NULL"}
    if(trim.blanks(powVal) != gettextRcmdr("")){
      nBlank <- nBlank + 1
      pow <- as.numeric(powVal)
      if(pow<=0 || pow>=1){
        errorCondition(recall=powerTtest, message=gettextRcmdr("Power must be between 0 and 1 if supplied."))
        return()
      }
    } else {powVal <- "NULL"}
    if(nBlank != 4){
      errorCondition(recall=powerTtest, message=gettextRcmdr("Exactly one field must be left empty."))
      return()
    }
    type <- as.character(tclvalue(typeVariable))
    alternative <- as.character(tclvalue(alternativeVariable))
    closeDialog()
    command <- paste("power.t.test(n = ", nVal, ", delta = ", deltaVal, ", sd = ", sdVal, ", sig.level = ", sigVal, ", power = ", powVal, ", type = '", type, "', alternative = '", alternative, "')", sep="")
    doItAndPrint(command)
    #assign(".Table", justDoIt(command), envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  nFrame <- tkframe(top)
  nLevel <- tclVar("")
  nField <- ttkentry(nFrame, width="6", textvariable=nLevel)
  deltaFrame <- tkframe(top)
  deltaLevel <- tclVar("")
  deltaField <- ttkentry(deltaFrame, width="6", textvariable=deltaLevel)
  sdFrame <- tkframe(top)
  sdLevel <- tclVar("")
  sdField <- ttkentry(sdFrame, width="6", textvariable=sdLevel)
  sigFrame <- tkframe(top)
  sigLevel <- tclVar("0.05")
  sigField <- ttkentry(sigFrame, width="6", textvariable=sigLevel)
  powFrame <- tkframe(top)
  powLevel <- tclVar("")
  powField <- ttkentry(powFrame, width="6", textvariable=powLevel)
  
  tkgrid(labelRcmdr(nFrame, text=gettextRcmdr("# of samples:"), fg="blue"), sticky="nw")
  tkgrid(nField, sticky="nw")
  tkgrid(labelRcmdr(deltaFrame, text=gettextRcmdr("True difference:"), fg="blue"), sticky="nw")
  tkgrid(deltaField, sticky="nw")
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("Standard deviation:"), fg="blue"), sticky="nw")
  tkgrid(sdField, sticky="nw")
  tkgrid(labelRcmdr(sigFrame, text=gettextRcmdr("Significance level:"), fg="blue"), sticky="nw")
  tkgrid(sigField, sticky="nw")
  tkgrid(labelRcmdr(powFrame, text=gettextRcmdr("Power:"), fg="blue"), sticky="nw")
  tkgrid(powField, sticky="nw")
  
  tkgrid(nFrame, sticky="nw", row=1, column=1, columnspan=1)
  tkgrid(deltaFrame, sticky="nw", row=1, column=2, columnspan=1)
  tkgrid(sdFrame, sticky="nw", row=1, column=3, columnspan=1)
  tkgrid(sigFrame, sticky="nw", row=2, column=1, columnspan=1)
  tkgrid(powFrame, sticky="nw", row=2, column=2, columnspan=1)
  
  radioButtons(top, name="type", buttons=c("twosample", "onesample", "Paired"), values=c("two.sample", "one.sample", "paired"),
               labels=gettextRcmdr(c("Two sample", "One sample", "Paired")), title=gettextRcmdr("Type of t test"))
  radioButtons(top, name="alternative", buttons=c("twosided", "onesided"), values=c("two.sided", "one.sided"),
               labels=gettextRcmdr(c("Two-sided", "One-sided")), title=gettextRcmdr("Alternative"))
  tkgrid(typeFrame, sticky="nw", row=3, column=1, columnspan=1)
  tkgrid(alternativeFrame, sticky="nw", row=3, column=2, columnspan=2)
  
  OKCancelHelp(helpSubject="power.t.test")
  tkgrid(buttonsFrame, sticky="w", row=4, column=1, columnspan=3)
  dialogSuffix(rows=4, columns=2)
}