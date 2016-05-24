fncCoinSpearmanTest <- function(){
  initializeDialog(title=gettextRcmdr("Spearman Test"))
  variablesFrame <- tkframe(top) #
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("First variable\n(select one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Second variable\n(select one)"))
  blockBox <- variableListBox(variablesFrame, Factors(), title="Block\n(select none or one)") #
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    subset <- tclvalue(subsetVariable) #
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
    test <- as.character(tclvalue(testVariable))
    if (length(x) == 0 | length(y) == 0) {
      errorCondition(recall=fncCoinSpearmanTest, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y) {
      errorCondition(recall=fncCoinSpearmanTest, message=gettextRcmdr("The two variables must be different."))
      return()
    }
    block <- getSelection(blockBox) #
    if (length(block) > 0) {
      block = paste(" | ", block, " ", sep="") #
    } else {
      block = "" #
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    
    if (test == "default"){
      doItAndPrint(paste("spearman_test(", x, " ~ ", y, block, ", data=", .activeDataSet, subset, ")", sep=""))
    }
    else doItAndPrint(paste("spearman_test(", x, " ~ ", y, block, ", distribution='", test, "'", ", data=", .activeDataSet, subset, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="spearman_test")
  optionsFrame <- tkframe(top) #
  
  radioButtons(optionsFrame, name="test", buttons=c("default", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test"))
  
  subsetBox() #
  
  tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox), getFrame(blockBox), sticky="nw") # 
  tkgrid(variablesFrame, sticky="nw") #
  tkgrid(testFrame, sticky="nw") #
  tkgrid(optionsFrame, sticky="nw") #
  #tkgrid(tiesFrame, sticky="w") #
  tkgrid(subsetFrame, sticky="w") #
  tkgrid(buttonsFrame, sticky="w") #
  dialogSuffix(rows=6, columns=1) #
}