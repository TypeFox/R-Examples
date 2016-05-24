
fncCoinPairedWilcoxonTest <- function(){
  initializeDialog(title=gettextRcmdr("Paired Wilcoxon Test"))
  variablesFrame <- tkframe(top) #
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("First variable\n(select one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Second variable\n(select one)"))
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    #exactoption <- as.character(tclvalue(exactbuttonsVariable))
    maxpts <- as.integer(tclvalue(asymptoticMaxpts))
    abseps <- as.double(tclvalue(asymptoticAbseps))
    releps <- as.double(tclvalue(asymptoticReleps))
    replications <- as.integer(tclvalue(approximateReplications))
    zeromethod <- as.character(tclvalue(zeromethodVariable))
    closeDialog()
    strExactOption = ""
    #if (exactoption == "default") {
    #    strExactOption = ""
    #} else {
    #    strExactOption = paste("algorithm='", exactoption, "'", sep="")
    #}
    strZeroMethod = ""
    if (zeromethod == "Default") {
      strZeroMethod = ""
    } else {
      strZeroMethod = paste(", zero.method='", zeromethod, "' ", sep="")
    }
    alternative <- as.character(tclvalue(alternativeVariable))
    subset <- tclvalue(subsetVariable) #
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
    test <- as.character(tclvalue(testVariable))
    if (length(x) == 0 | length(y) == 0) {
      errorCondition(recall=fncCoinPairedWilcoxonTest, message=gettextRcmdr("You must select two variables."))
      return()
    }
    if (x == y) {
      errorCondition(recall=fncCoinPairedWilcoxonTest, message=gettextRcmdr("The two variables must be different."))
      return()
    }
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("median(", .activeDataSet, "$", x, " - ", .activeDataSet, "$", y, 
                       ", na.rm=TRUE) # median difference", sep=""))
    if (test == "default"){
      doItAndPrint(paste("wilcoxsign_test(", .activeDataSet, "$", x, " ~ ", .activeDataSet, "$", y,
                         ", alternative='",  alternative, "' ", strZeroMethod, subset, ")", sep=""))           
    } else if (test == "exact"){
      doItAndPrint(paste("wilcoxsign_test(", .activeDataSet, "$", x, " ~ ", .activeDataSet, "$", y,
                         ", alternative='", alternative, "' ", strZeroMethod, subset, ", distribution=exact(", strExactOption, "))", sep=""))
    } else if (test == "asympt"){
      doItAndPrint(paste("wilcoxsign_test(", .activeDataSet, "$", x, " ~ ", .activeDataSet, "$", y,
                         ", alternative='", alternative, "' ", strZeroMethod, subset, ", distribution=asymptotic(maxpts = ", maxpts, ", abseps = ", abseps, ", releps = ", releps, " ) )", sep=""))
    } else {
      doItAndPrint(paste("wilcoxsign_test(", .activeDataSet, "$", x, " ~ ", .activeDataSet, "$", y,
                         ", alternative='", alternative, "' ", strZeroMethod, subset, ", distribution=approximate(B = ", replications, " ) )", sep=""))
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="wilcox_test")
  optionsFrame <- tkframe(top) #
  radioButtons(optionsFrame, name="alternative", buttons=c("two.sided", "less", "greater"), 
               labels=gettextRcmdr(c("Two-sided", "Difference < 0", "Difference > 0")), title=gettextRcmdr("Alternative Hypothesis"))
  radioButtons(optionsFrame, name="test", buttons=c("default", "exact", "asympt", "approx"), 
               labels=c("Default", "Exact", "Asymptotic", "Approximate"), 
               title=gettextRcmdr("Type of Test"))
  
  distributionsFrame <- tkframe(top)
  #radioButtons(distributionsFrame, name="exactbuttons", buttons=c("Default", "Shift", "SplitUp"), values=c("default", "shift", "split-up"), labels=c("Default", "Shift", "Split-up"), title="Exact options")
  
  asymptoticFrame <- tkframe(distributionsFrame) #
  asymptoticMaxptsFrame <- tkframe(asymptoticFrame) #
  asymptoticMaxpts <- tclVar("25000") #
  asymptoticMaxptsField <- ttkentry(asymptoticMaxptsFrame, width="6", textvariable=asymptoticMaxpts) #
  asymptoticAbsepsFrame <- tkframe(asymptoticFrame) #
  asymptoticAbseps <- tclVar("0.001") #
  asymptoticAbsepsField <- ttkentry(asymptoticAbsepsFrame, width="6", textvariable=asymptoticAbseps) #
  asymptoticRelepsFrame <- tkframe(asymptoticFrame) #
  asymptoticReleps <- tclVar("0") #
  asymptoticRelepsField <- ttkentry(asymptoticRelepsFrame, width="6", textvariable=asymptoticReleps) #
  
  approximateFrame <- tkframe(distributionsFrame) #
  approximateReplications <- tclVar("1000") #
  approximateField <- ttkentry(approximateFrame, width="6", textvariable=approximateReplications) #
  
  zeroMethodFrame <- tkframe(top) #
  radioButtons(zeroMethodFrame, name="zeromethod", buttons=c("Default", "Pratt", "Wilcoxon"), 
               labels=c("Default", "Pratt", "Wilcoxon"), title="Zero method")
  
  subsetBox() #
  
  tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw") #
  tkgrid(variablesFrame, sticky="nw") #
  
  tkgrid(alternativeFrame, labelRcmdr(optionsFrame, text="    "), testFrame, sticky="nw") #
  tkgrid(optionsFrame, sticky="nw") #
  
  tkgrid(labelRcmdr(distributionsFrame, text="Distribution options", fg="blue"),sticky="w") #
  
  
  #tkgrid(exactbuttonsFrame, sticky="w") #
  
  
  tkgrid(labelRcmdr(asymptoticFrame, text="Asymptotic options", fg="blue"),sticky="w") #
  tkgrid(labelRcmdr(asymptoticMaxptsFrame, text="Max pts", fg="blue"),sticky="w") #
  tkgrid(asymptoticMaxptsField, sticky="w") #
  tkgrid(labelRcmdr(asymptoticAbsepsFrame, text="Abs. err. tolerance", fg="blue"),sticky="w") #
  tkgrid(asymptoticAbsepsField, sticky="w") #
  tkgrid(labelRcmdr(asymptoticRelepsFrame, text="Rel. err. tolerance", fg="blue"),sticky="w") #
  tkgrid(asymptoticRelepsField, sticky="w") #
  
  tkgrid(asymptoticMaxptsFrame, asymptoticAbsepsFrame, asymptoticRelepsFrame, sticky="s") #
  tkgrid(asymptoticFrame, sticky="nw")
  
  tkgrid(labelRcmdr(approximateFrame, text="Approximate replications", fg="blue"),sticky="w") #
  tkgrid(approximateField, sticky="w") #
  
  tkgrid(approximateFrame, labelRcmdr(distributionsFrame, text="    "), sticky="nw") #
  tkgrid(distributionsFrame, sticky="nw") #
  
  tkgrid(zeromethodFrame, labelRcmdr(zeroMethodFrame, text="    "), sticky="nw") #
  tkgrid(zeroMethodFrame, sticky="nw") #
  
  tkgrid(subsetFrame, sticky="w") #
  tkgrid(buttonsFrame, sticky="w") #
  dialogSuffix(rows=6, columns=1) #
  #    tkgrid(getFrame(xBox), getFrame(yBox), sticky="nw")    
  #    tkgrid(alternativeFrame, testFrame, sticky="nw")
  #    tkgrid(buttonsFrame, columnspan=2, sticky="w")
  #    dialogSuffix(rows=3, columns=1)
}