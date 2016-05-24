# Statistics Menu dialogs

# last modified 2014-08-17 by J. Fox

# Dimensional-analysis menu

Reliability <- function () {
	defaults <- list(initial.x = NULL)
	dialog.values <- getDialog("Reliability", defaults)
	initializeDialog(title = gettextRcmdr("Scale Reliability"))
	xBox <- variableListBox(top, Numeric(), selectmode = "multiple", 
			initialSelection = varPosn(dialog.values$initial.x, "numeric"),
			title = gettextRcmdr("Variables (pick three or more)"))
	onOK <- function() {
		x <- getSelection(xBox)
		closeDialog()
		putDialog("Reliability", list (initial.x = x))
		if (3 > length(x)) {
			errorCondition(recall = Reliability, message = gettextRcmdr("Fewer than 3 variables selected."))
			return()
		}
		x <- paste("\"", x, "\"", sep = "")
		doItAndPrint(paste("reliability(cov(", ActiveDataSet(), 
						"[,c(", paste(x, collapse = ","), ")], use=\"complete.obs\"))", 
						sep = ""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "reliability", reset = "Reliability", apply = "Reliability")
	tkgrid(getFrame(xBox), sticky = "nw")
	tkgrid(buttonsFrame, sticky = "w")
	dialogSuffix()
}

principalComponents <- function () {
  defaults <- list(initial.x = NULL, initial.correlations = 1, 
                   initial.subset = gettextRcmdr("<all valid cases>"), initial.screeplot = 0, initial.addPC = 0,
                   initial.tab=0)
  dialog.values <- getDialog("principalComponents", defaults)
  initializeDialog(title = gettextRcmdr("Principal Components Analysis"), use.tabs=TRUE)
  xBox <- variableListBox(dataTab, Numeric(), selectmode = "multiple", 
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"), 
                          title = gettextRcmdr("Variables (pick two or more)"))
  subsetBox(dataTab, subset.expression = dialog.values$initial.subset)
  checkBoxes(optionsTab, frame = "optionsFrame", boxes = c("correlations", 
                                                           "screeplot", "addPC"), initialValues = c(dialog.values$initial.correlations, 
                                                                                                    dialog.values$initial.screeplot, dialog.values$initial.addPC), 
             labels = gettextRcmdr(c("Analyze correlation matrix", 
                                     "Screeplot", "Add principal components to data set")))
  onOK <- function() {
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    putRcmdr("ncomponents", 0)
    x <- getSelection(xBox)
    nvar <- length(x)
    correlations <- tclvalue(correlationsVariable)
    subset <- tclvalue(subsetVariable)
    screeplot <- tclvalue(screeplotVariable)
    addPC <- tclvalue(addPCVariable)
    closeDialog()
    putDialog("principalComponents", list(initial.x = x, initial.correlations = correlations, 
                                          initial.subset = subset, initial.screeplot = screeplot, initial.addPC = addPC,
                                          initial.tab=tab))
    if (2 > length(x)) {
      errorCondition(recall = principalComponents, message = gettextRcmdr("Fewer than 2 variables selected."))
      return()
    }
    subset <- if (trim.blanks(subset) == "" || trim.blanks(subset) == gettextRcmdr("<all valid cases>")) 
      ""
    else paste(", subset=", subset, sep = "")
    correlations <- if (correlations == "1") 
      "TRUE"
    else "FALSE"
    .activeDataSet <- ActiveDataSet()
    cmds <- character(8)
    cmds[1] <- paste("local({\n  .PC <- princomp(~", paste(x, collapse = "+"), 
                     ", cor=", correlations, ", data=", .activeDataSet, 
                     subset, ")", sep = "")
    cmds[2] <- '  cat("\\nComponent loadings:\\n")'
    cmds[3] <- "  print(unclass(loadings(.PC)))"
    cmds[4] <- '  cat("\\nComponent variances:\\n")'
    cmds[5] <- "  print(.PC$sd^2)"
    cmds[6] <- '  cat("\\n")'
    cmds[7] <- "  print(summary(.PC))"
    cmds[8] <- if (screeplot == "1") "  screeplot(.PC)\n" else ""
    cmds <- paste(cmds[1:(7 + (screeplot == "1"))], collapse="\n")
    commands <- ""
    if (addPC == "1") {
      if (trim.blanks(subset) != ""){
        errorCondition(recall=principalComponents,
                       message=gettextRcmdr("Component scores are not available when subset is specified."))
        return()
      }
      initializeDialog(subdialog, title = gettextRcmdr("Number of Components"))
      tkgrid(labelRcmdr(subdialog, text = gettextRcmdr("Number of components to retain:"), 
                        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
      sliderFrame <- tkframe(subdialog)
      sliderValue <- tclVar("1")
      componentsSlider <- tkscale(sliderFrame, from = 1, 
                                  to = nvar, showvalue = FALSE, variable = sliderValue, 
                                  resolution = 1, orient = "horizontal")
      componentsShow <- labelRcmdr(sliderFrame, textvariable = sliderValue, 
                                   width = 2, justify = "right")
      onOKsub <- function() {
        closeDialog(subdialog)
        putRcmdr("ncomponents", as.numeric(tclvalue(sliderValue)))
      }
      subOKCancelHelp()
      tkgrid(componentsSlider, componentsShow, sticky = "nw")
      tkgrid(sliderFrame, sticky = "w")
      tkgrid(subButtonsFrame, sticky = "w")
      dialogSuffix(subdialog, onOK = onOKsub, focus = subdialog, force.wait=TRUE)
      if ((ncomponents <- getRcmdr("ncomponents")) > 0) {
        prefix <- paste("  ", .activeDataSet, " <<- within(", .activeDataSet, ", {", sep="")
        if (screeplot != "1") prefix <- paste("\n", prefix, sep="")
        commands <- character(ncomponents)
        for (i in 1:ncomponents) {
          var <- paste("PC", i, sep = "")
          if (is.element(var, Variables())) {
            if ("no" == tclvalue(checkReplace(var))) 
              next
          }
          commands[ncomponents - i + 1] <- paste("    PC", i, " <- .PC$scores[,", 
                                                 i, "]", sep = "")
        }
        suffix <- "  })"
        commands <- paste(c(prefix, commands, suffix), collapse="\n")
      }
    }
    doItAndPrint(paste(cmds, commands, "\n})", sep=""))
    if (addPC == "1") activeDataSet(.activeDataSet, flushDialogMemory=FALSE)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "princomp", reset = "principalComponents", apply = "principalComponents")
  tkgrid(getFrame(xBox), sticky = "nw")
  tkgrid(subsetFrame, sticky = "w")
  tkgrid(optionsFrame, sticky = "w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

factorAnalysis <- function () {
  defaults <- list(initial.x = NULL, initial.subset = gettextRcmdr ("<all valid cases>"), 
                   initial.rotation = "varimax", initial.scores = "none", initial.tab=0)
  dialog.values <- getDialog("factorAnalysis", defaults)
  initializeDialog(title = gettextRcmdr("Factor Analysis"), use.tabs=TRUE)
  xBox <- variableListBox(dataTab, Numeric(), selectmode = "multiple", 
                          initialSelection = varPosn(dialog.values$initial.x, "numeric"),
                          title = gettextRcmdr("Variables (pick three or more)"))
  subsetBox(dataTab, subset.expression = dialog.values$initial.subset)
  optionsFrame <- tkframe(optionsTab)
  checkFrame <- tkframe(optionsTab)
  radioButtons(checkFrame, name = "rotation", buttons = c("noRotate", 
                                                          "varimax", "promax"), values = c("none", "varimax", "promax"), 
               initialValue = dialog.values$initial.rotation, labels = gettextRcmdr(c("None", 
                                                                                      "Varimax", "Promax")), title = gettextRcmdr("Factor Rotation"))
  radioButtons(checkFrame, name = "scores", buttons = c("noScores", 
                                                        "bartlett", "regression"), values = c("none", "Bartlett", 
                                                                                              "regression"), initialValue = dialog.values$initial.scores,
               labels = gettextRcmdr(c("None", "Bartlett's method", 
                                       "Regression method")), title = gettextRcmdr("Factor Scores"))
  onOK <- function() {
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(xBox)
    nvar <- length(x)
    subset <- tclvalue(subsetVariable)
    rotation <- tclvalue(rotationVariable)
    scores <- tclvalue(scoresVariable)
    closeDialog()
    putDialog ("factorAnalysis", list (initial.x = x, initial.subset = subset, 
                                       initial.scores = scores, initial.rotation = rotation, initial.tab=tab))
    if (3 > length(x)) {
      errorCondition(recall = factorAnalysis, message = gettextRcmdr("Fewer than 3 variables selected."))
      return()
    }
    f <- function(k, p) ((p - k)^2 - p - k)^2
    max.factors <- floor(optimize(f, c(0, nvar), tol = 1e-04, 
                                  p = nvar)$minimum)
    if (max.factors == 1) {
      putRcmdr("nfactors", 1)
    }
    else {
      initializeDialog(subdialog, title = gettextRcmdr("Number of Factors"))
      tkgrid(labelRcmdr(subdialog, text = gettextRcmdr("Number of factors to extract:"), 
                        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
      sliderFrame <- tkframe(subdialog)
      sliderValue <- tclVar("1")
      componentsSlider <- tkscale(sliderFrame, from = 1, 
                                  to = max.factors, showvalue = FALSE, variable = sliderValue, 
                                  resolution = 1, orient = "horizontal")
      componentsShow <- labelRcmdr(sliderFrame, textvariable = sliderValue, 
                                   width = 2, justify = "right")
      onOKsub <- function() {
        closeDialog(subdialog)
        putRcmdr("nfactors", as.numeric(tclvalue(sliderValue)))
      }
      subOKCancelHelp()
      tkgrid(componentsSlider, componentsShow, sticky = "nw")
      tkgrid(sliderFrame, sticky = "w")
      tkgrid(subButtonsFrame, sticky = "w")
      dialogSuffix(subdialog, onOK = onOKsub, focus = subdialog, force.wait=TRUE)
    }
    subset <- if (trim.blanks(subset) == "" || trim.blanks(subset) == gettextRcmdr("<all valid cases>")) 
      ""
    else paste(", subset=", subset, sep = "")
    if (scores != "none" && subset != ""){
      errorCondition(recall=factorAnalysis,
                     message=gettextRcmdr("Factor scores are not available when subset is specified."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    command <- paste("local({\n  .FA <- factanal(~", paste(x, collapse = "+"), 
                     ", factors=", getRcmdr("nfactors"), ", rotation=\"", 
                     rotation, "\", scores=\"", scores, "\", data=", .activeDataSet, 
                     subset, ")\n  print(.FA)\n", sep = "")
    commands <- ""
    if (scores != "none") {
      prefix <- paste("  ", .activeDataSet, " <<- within(", .activeDataSet, ", {", sep="")
      nfactors <- getRcmdr("nfactors")
      commands <- character(nfactors)
      for (i in 1:nfactors) {
        var <- paste("F", i, sep = "")
        if (is.element(var, Variables())) {
          if ("no" == tclvalue(checkReplace(var))) 
            next
        }
        commands[nfactors - i + 1] <- paste("    ", var, " <- .FA$scores[,", i, "]", sep = "")
      }
      suffix <- "  })"
      commands <- paste(c(prefix, commands, suffix), collapse="\n")
    }
    if (commands != "") doItAndPrint(paste(command, commands, "\n})", sep=""))
    else doItAndPrint(paste(command, "})", sep=""))
    if (scores != "none") activeDataSet(.activeDataSet, flushDialogMemory = FALSE)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "factanal", reset = "factorAnalysis", apply = "factorAnalysis")
  tkgrid(getFrame(xBox), sticky = "nw")
  tkgrid(subsetFrame, sticky = "w")
  tkgrid(optionsFrame, sticky = "w")
  tkgrid(rotationFrame, labelRcmdr(checkFrame, text = "    "), 
         scoresFrame, sticky = "w")
  tkgrid(checkFrame, sticky = "w")
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

CFA <- function(){
  Library("sem")
  defaults <- list(initial.matrix="covariance", initial.factorCor="correlated", initial.identify="factors", 
                   initial.robust=0, initial.tab=0,
                   initial.AIC=1, initial.BIC=1, initial.GFI=0, initial.AGFI=0,
                   initial.RMSEA=0, initial.NFI=0, initial.NNFI=0,
                   initial.CFI=0, initial.RNI=0, initial.IFI=0, initial.SRMR=0,
                   initial.AICc=0, initial.CAIC=0)
  dialog.values <- getDialog("CFA", defaults)
  initializeDialog(title=gettextRcmdr("Confirmatory Factor Analysis"), use.tabs=TRUE)
  onFactor <- function(){
    vars <- getSelection(xBox)
    if (length(vars) < 2) {
      errorCondition(recall=CFA,  message=gettextRcmdr("Fewer than 2 variables selected to load on factor."))
      return()
    }
    fac.name <- tclvalue(factorName)
    if (!is.valid.name(fac.name)) {
      errorCondition(recall=CFA,  message=paste(fac.name, gettextRcmdr("is not a valid name.")))
      return()
    }
    variables[[getRcmdr("factorNumber")]] <<- vars
    factors <- factors[getRcmdr("factorNumber")] <<- fac.name
    putRcmdr("factorNumber", getRcmdr("factorNumber") + 1)
    tclvalue(factorName) <- paste("Factor.", getRcmdr("factorNumber"), sep = "")
    tkselection.clear(xBox$listbox, "0", "end")
    tclvalue(buttonText) <- paste(gettextRcmdr("Define factor"), getRcmdr("factorNumber"))
  }
  xBox <- variableListBox(dataTab, Numeric(), selectmode = "multiple", 
                          title = gettextRcmdr("Select variables\nloading on factor"))
  optionsFrame <- tkframe(optionsTab)
  radioButtons(optionsFrame, name = "matrix", buttons = c("covariance", "correlation"),
               initialValue = dialog.values$initial.matrix, 
               labels = gettextRcmdr(c("Covariance", "Correlation")), title = gettextRcmdr("Matrix to Analyze"))
  radioButtons(optionsFrame, name = "factorCor", buttons = c("correlated", "orthogonal"),
               initialValue = dialog.values$initial.factorCor, 
               labels = gettextRcmdr(c("Correlated", "Orthogonal")), title = gettextRcmdr("Factor Correlations"))
  radioButtons(optionsFrame, name = "identify", buttons = c("factors", "loadings"),
               initialValue = dialog.values$initial.identify, 
               labels = gettextRcmdr(c("Factor variances set to 1", "First loading on each factor set to 1")), 
               title = gettextRcmdr("Identifying Constraints"))
  fit.indices.1 <- c("AIC", "BIC", "GFI", "AGFI", "RMSEA", "NFI", "NNFI")
  fit.indices.2 <- c("CFI", "RNI", "IFI", "SRMR", "AICc", "CAIC")
  fit.indices <- c(fit.indices.1, fit.indices.2)
  fitIndicesFrame <- tkframe(optionsTab)
  checkBoxes(window=fitIndicesFrame, frame="fitIndicesFrame.1", 
             boxes=fit.indices.1,
             initialValues=unlist(dialog.values[paste("initial.", fit.indices.1, sep="")]),
             labels=fit.indices.1,
             title=gettextRcmdr("Fit Indices"))
  checkBoxes(window=fitIndicesFrame, frame="fitIndicesFrame.2", 
             boxes=fit.indices.2,
             initialValues=unlist(dialog.values[paste("initial.", fit.indices.2, sep="")]),
             labels=fit.indices.2,
             title="")
  checkBoxes(window=optionsFrame, frame = "robustFrame", boxes = "robust", initialValues = dialog.values$initial.robust, 
             labels = gettextRcmdr("Robust standard errors"), title=" ")
  putRcmdr("factorNumber", 1)
  buttonText <- tclVar(paste(gettextRcmdr("Define factor"), getRcmdr("factorNumber")))
  factorFrame <- tkframe(dataTab)
  factorButton <- buttonRcmdr(factorFrame, textvariable=buttonText, width="15", 
                              command=onFactor, default="active", borderwidth=3)
  factorName <- tclVar(paste("Factor.", getRcmdr("factorNumber"), sep = ""))
  factorEntry <- ttkentry(factorFrame, width="20", textvariable=factorName)
  variables <- list()
  factors <- vector()
  onOK <- function(){
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    matrix <- tclvalue(matrixVariable)
    correlations <- tclvalue(factorCorVariable)
    identify <- tclvalue(identifyVariable)
    robust <- tclvalue(robustVariable)
    for (index in fit.indices){
      assign(index, tclvalue(eval(parse(text=paste(index, "Variable", sep="")))))
    }
    indices <- fit.indices[as.logical(as.numeric(c(AIC, BIC, GFI, AGFI, RMSEA, NFI, 
                                                   NNFI, CFI, RNI, IFI, SRMR, AICc, CAIC)))]
    indices <- paste(', fit.indices=c(', paste(paste('"', indices, '"', sep=""), collapse=","),
                     ')', sep="")
    closeDialog()
    putDialog("CFA", list(initial.matrix=matrix, initial.factorCor=correlations, 
                          initial.identify=identify, initial.robust=robust, initial.tab=tab,
                          initial.AIC=AIC, initial.BIC=BIC, initial.GFI=GFI, initial.AGFI=AGFI,
                          initial.RMSEA=RMSEA, initial.NFI=NFI, initial.NNFI=NNFI,
                          initial.CFI=CFI, initial.RNI=RNI, initial.IFI=IFI, initial.SRMR=SRMR,
                          initial.AICc=AICc, initial.CAIC=CAIC))
    if (length(factors) == 0) {
      errorCondition(recall=CFA, message=gettextRcmdr("No factors defined."))
      return()
    }
    putRcmdr("factorNumber", NULL)
    modelText <- vector(length(factors), mode="character")
    for (i in 1:length(factors)){
      modelText[i] <- paste(factors[i], ": ", paste(variables[[i]], collapse=", "), sep="")
    }
    allvars <- unique(unlist(variables))
    if ((length(allvars)/length(factors)) < 2) {
      errorCondition(recall=CFA,  
                     message=gettextRcmdr("There are too many factors."))
      return()
    }
    command <- paste("local({\n  .model <- c(", paste(paste("'", modelText, "'", sep=""), 
                                                      collapse=", "), ")", sep="")
    command <- paste(command, "\n  .model <- cfa(file=textConnection(.model)", 
                     if(correlations == "correlated") ", " else ", covs=NULL, ",
                     "reference.indicators=", 
                     if (identify == "factors") "FALSE" else "TRUE", ")", sep="")
    command <- paste(command, "\n  .Data <- ", activeDataSet(),
                     "[, c(", paste(paste("'", allvars, "'", sep=""), collapse=", "),  ")]", 
                     sep="")
    if (matrix == "correlation") 
      command <- paste(command, "\n  .Data <- as.data.frame(scale(.Data))")
    command <- paste(command, "\n  summary(sem(.model, data=.Data), robust=", 
                     if (robust == 1) "TRUE" else "FALSE", indices,
                     ")\n})", sep="")
    doItAndPrint(command)
  }
  OKCancelHelp(helpSubject="CFA", reset="CFA", apply="CFA")
  tkgrid(matrixFrame, labelRcmdr(optionsFrame, text="    "), factorCorFrame, sticky="nw")
  tkgrid(identifyFrame, labelRcmdr(optionsFrame, text="    "),  robustFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(fitIndicesFrame.1, fitIndicesFrame.2, sticky="nw")
  tkgrid(fitIndicesFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="w")
  tkgrid(labelRcmdr(dataTab, text=""))
  tkgrid(factorButton, labelRcmdr(factorFrame, text=paste("   ", gettextRcmdr("Name for factor:"))), factorEntry, sticky="nw")
  tkgrid(factorFrame, sticky="w")
  tkgrid(labelRcmdr(dataTab, text=""))
  dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

