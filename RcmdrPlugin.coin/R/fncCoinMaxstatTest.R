fncCoinMaxstatTest <- function(){
  initializeDialog(title=gettextRcmdr("Maximally Selected Statistics Test"))
  variablesFrame <- tkframe(top) #
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, selectmode="multiple",
                          title=gettextRcmdr("Explanatory variables\n(select one or more)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Response variable\n(select one)"))
  blockBox <- variableListBox(variablesFrame, Factors(), title="Block\n(select none or one)") #
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    
    if (0 == length(y)) {
      errorCondition(recall=fncCoinMaxstatTest, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    if (0 == length(x)) {
      errorCondition(recall=fncCoinMaxstatTest, message=gettextRcmdr("No explanatory variables selected."))
      return()
    }
    if (is.element(y, x)) {
      errorCondition(recall=fncCoinMaxstatTest, message=gettextRcmdr("Response and explanatory variables must be different."))
      return()
    }
    block <- getSelection(blockBox) #
    if (length(block) > 0) {
      block = paste(" | ", block, " ", sep="") #
    } else {
      block = "" #
    }
    
    test <- as.character(tclvalue(testVariable)) #
    subset <- tclvalue(subsetVariable) #
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    if (test == "default"){
      doItAndPrint(paste("maxstat_test(", y, " ~ ", paste(x, collapse="+"), block, ", data=", .activeDataSet, subset, ")", sep=""))
    }
    else doItAndPrint(paste("maxstat_test(", y, " ~ ", paste(x, collapse="+"), block, 
                            ", distribution='", test, "'", ", data=", .activeDataSet, subset, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="maxstat_test")
  optionsFrame <- tkframe(top) #
  
  radioButtons(optionsFrame, name="test", buttons=c("default", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test"))
  
  subsetBox() #
  
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), getFrame(blockBox), sticky="nw") # 
  tkgrid(variablesFrame, sticky="nw") #
  tkgrid(testFrame, sticky="nw") #
  tkgrid(optionsFrame, sticky="nw") #
  tkgrid(subsetFrame, sticky="w") #
  tkgrid(buttonsFrame, sticky="w") #
  dialogSuffix(rows=6, columns=1) #
}