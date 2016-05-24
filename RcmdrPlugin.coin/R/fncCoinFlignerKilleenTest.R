fncCoinFlignerKilleenTest <- function(){
  initializeDialog(title=gettextRcmdr("Fligner Killeen Test"))
  variablesFrame <- tkframe(top) #
  groupBox <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Groups\n(select one)"))
  responseBox <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Response Variable\n(select one)"))
  blockBox <- variableListBox(variablesFrame, Factors(), title="Block\n(select none or one)") #
  onOK <- function(){
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall=fncCoinFlignerKilleenTest , message=gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall=fncCoinFlignerKilleenTest , message=gettextRcmdr("You must select a response variable."))
      return()
    }
    
    #level <- tclvalue(confidenceLevel) #
    #strConfintText <- paste (", conf.int = TRUE, conf.level=", level , sep="") #
    #strConfintText <- "conf.int = FALSE"
    #vties <- as.character(tclvalue(tiesVariable))
    
    block <- getSelection(blockBox) #
    if (length(block) > 0) {
      if (block == group) {
        errorCondition(recall=fncCoinFlignerKilleenTest, message=gettextRcmdr("The group and block variables must be different."))
        return()
      }
      block = paste(" | ", block, " ", sep="") #
      strConfintText = "" # cannot compute test for blocks !!!
    } else {
      block = "" #
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    
    test <- as.character(tclvalue(testVariable)) #
    subset <- tclvalue(subsetVariable) #
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    Library("abind")
    doItAndPrint(paste("numSummary(", paste(.activeDataSet,"$", response, sep=""),
                       ", ", paste("groups=",.activeDataSet,"$", group, sep=""), ", statistics=c('quantiles', 'sd'), quantiles=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))", sep=""))
    Ties <- as.character(tclvalue(tiesmethodVariable)) #
    if (Ties == "default") {
      strTies = ""
    } else {
      strTies = paste(", ties.method='", Ties, "' ", sep="")
    }
    if (test == "default"){
      doItAndPrint(paste("fligner_test(", response, " ~ ", group, block, ", alternative='", 
                         alternative, "'", strTies, ", data=", .activeDataSet, subset, ")", sep=""))
    }
    else doItAndPrint(paste("fligner_test(", response, " ~ ", group, block, ", alternative='", 
                            alternative, "'", ", distribution='", test, "'", strTies, ", data=", .activeDataSet, subset, ")", sep=""))
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="fligner_test")
  optionsFrame <- tkframe(top) #
  radioButtons(optionsFrame, name="alternative", buttons=c("twosided", "less", "greater"), values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Two-sided", "Difference < 0", "Difference > 0")), title=gettextRcmdr("Alternative Hypothesis"))
  
  #confidenceFrame <- tkframe(optionsFrame) #
  #confidenceLevel <- tclVar(".95") #
  #confidenceField <- ttkentry(confidenceFrame, width="6", textvariable=confidenceLevel) #
  
  radioButtons(optionsFrame, name="test", buttons=c("default", "approximate", "asymptotic"), 
               labels=gettextRcmdr(c("Default", "Monte Carlo resampling approximation", "Asymptotic null distribution")), 
               title=gettextRcmdr("Type of Test"))
  
  radioButtons(optionsFrame, name="tiesmethod",
               buttons=c("default", "midranks", "averagescores"), #
               values=c("default", "mid-ranks", "average-scores"), initialValue="default",
               labels=c("Default", "Mid-ranks", "Average scores"), 
               title=gettext("Method for ties"))	    
  
  subsetBox() #
  
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(responseBox), getFrame(blockBox), sticky="nw") # getFrame(blockBox), 
  tkgrid(variablesFrame, sticky="nw") #
  #tkgrid(labelRcmdr(confidenceFrame, text=gettextRcmdr("Confidence Level"), fg="blue"),sticky="w") #
  #tkgrid(confidenceField, sticky="w") #
  groupsLabel(groupsBox=groupBox) #
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "),  labelRcmdr(optionsFrame, text="    "), testFrame, sticky="nw") #
  tkgrid(optionsFrame, sticky="nw") #
  tkgrid(tiesmethodFrame, sticky="new")
  tkgrid(subsetFrame, sticky="w") #
  tkgrid(buttonsFrame, sticky="w") #
  dialogSuffix(rows=6, columns=1) #
}