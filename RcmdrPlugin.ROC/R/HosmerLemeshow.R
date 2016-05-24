fndHosmerLemeshow <- function(){
  defaults <- list(initial.g = "10", 
                   initial.classifTable = 1, initial.percentCorrect = 1,
                   initial.tab=0)
  dialog.values <- getDialog("HosmerLemeshow", defaults)
  initializeDialog(title = gettext("Hosmer-Lemeshow GOF Test", domain="R-RcmdrPlugin.ROC"), use.tabs=TRUE)
  
  gVar <- tclVar(dialog.values$initial.g)
  gEntry <- ttkentry(dataTab, width = "25", textvariable = gVar)
  
  optionsFrame <- tkframe(optionsTab)# tab
  checkBoxes(window = optionsFrame, 
             frame = "aditionaloptionsFrame", # tab
             boxes = c("classifTable", "percentCorrect"), 
             initialValues = c(dialog.values$initial.classifTable, dialog.values$initial.percentCorrect),
             labels = gettextRcmdr(c("Classification table", "Percentage correct")), 
             title = gettextRcmdr(""), ttk=TRUE)
  
  onOK <- function() {
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    g <- as.numeric(as.character(tclvalue(gVar)))
    
    classifTable <- as.character("1" == tclvalue(classifTableVariable))
    percentCorrect <- as.character("1" == tclvalue(percentCorrectVariable))
    
    if (length(g) == 0) {
      errorCondition(recall = fndHosmerLemeshow, message = gettext("You must enter the number of bins to use to calculate quantiles.", domain="R-RcmdrPlugin.ROC"))
      return()
    }

    closeDialog()
    putDialog("HosmerLemeshow", list(initial.g = g, 
                                     initial.classifTable = tclvalue(classifTableVariable), initial.percentCorrect = tclvalue(percentCorrectVariable),
                                            initial.tab=tab))
    .activeDataSet <- ActiveDataSet()

    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) return()
     
    command <- paste(".depname <- ", "as.character((attr(", .activeModel, "$terms, 'variables')[2]))", sep = "")
    justDoIt(command)
    command <- paste(".outcome <- ifelse(", .activeDataSet, "$", .depname, "==levels(as.factor(", .activeDataSet, "$", .depname, "))[2], 1, 0 )", sep = "")
    doItAndPrint(command)
    # remove incomplete cases
    command <- paste(".matrix <- cbind (.outcome, fitted(", .activeModel, "))", sep = "")
    doItAndPrint(command)
    command <- paste(".matrix <- .matrix[complete.cases(.matrix),]", sep = "")
    doItAndPrint(command)
    command <- paste(".hl <-hoslem.test(.matrix[,1], .matrix[,2], ", g, ")", sep = "")
    doItAndPrint(command)
    command <- paste(".hl", sep = "")
    doItAndPrint(command)

    
    command <- paste(".treshold  <- 0.5", sep = "")
    doItAndPrint(command)
    command <- paste(".predictedBinary <- cut(fitted(", .activeModel, "), breaks=c(-Inf, .treshold, Inf), labels=c('low', 'high'))", sep = "")
    doItAndPrint(command)
    
    command <- paste(".tableClassif <- table(.outcome, .predictedBinary) # Classification table:", sep = "")
    doItAndPrint(command)
    if (classifTable == "TRUE") {
      command <- paste(".tableClassif", sep = "")
      doItAndPrint(command)
    }
    if (percentCorrect == "TRUE") {
      command <- paste(".percentageCorrect <- round(sum(diag(.tableClassif)) * 100 / sum(.tableClassif), 2) # Overall percentage correct:", sep = "")
      doItAndPrint(command)
      command <- paste(".percentageCorrect", sep = "")
      doItAndPrint(command)
    }
    
    # remove variables
    command <- paste("remove(.depname)", sep = "")
    justDoIt(command)
    command <- paste("remove(.outcome)", sep = "")
    doItAndPrint(command)
    command <- paste("remove(.matrix)", sep = "")
    doItAndPrint(command)
    command <- paste("remove(.hl)", sep = "")
    doItAndPrint(command)
    if (classifTable == "TRUE") {
      command <- paste("remove(.tableClassif)", sep = "")
      doItAndPrint(command)
    }
    if (percentCorrect == "TRUE") {
      command <- paste("remove(.percentageCorrect)", sep = "")
      doItAndPrint(command)
    }
    
#    command <- paste(".dep <- ", .activeDataSet, "$", "as.character((attr(", .activeModel, "$terms, 'variables')[2]))", sep = "")
    
    
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "hoslem.test", reset = "fndHosmerLemeshow",
               apply = "fndHosmerLemeshow")

tkgrid(aditionaloptionsFrame, sticky = "we")
tkgrid(optionsFrame, sticky = "we")

  tkgrid(labelRcmdr(dataTab, text = gettext("Number of bins", domain="R-RcmdrPlugin.ROC")), gEntry, sticky = "ew", padx=6, pady=c(6, 0))


  #dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("dataTab", "optionsTab"), 
             tab.names=c("Test options", "Classification")) #
}