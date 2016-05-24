fncCoinFriedmanTest <- function(){
  #Library("multcomp")
  initializeDialog(title=gettextRcmdr("Friedman rank-sum Test"))
  variablesFrame <- tkframe(top) #
  groupBox <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Groups\n(select one)"))
  responseBox <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Response Variable\n(select one)"))
  blockBox <- variableListBox(variablesFrame, Factors(), title="Block\n(select none or one)") #
  
  onOK <- function(){
    group <- getSelection(groupBox)
    if (length(group) == 0) {
      errorCondition(recall=fncCoinFriedmanTest, message=gettextRcmdr("You must select a groups variable."))
      return()
    }
    response <- getSelection(responseBox)
    if (length(response) == 0) {
      errorCondition(recall=fncCoinFriedmanTest, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    strBlock <- getSelection(blockBox) #
    if (length(strBlock) > 0) {
      if (strBlock == group) {
        errorCondition(recall=fncCoinFriedmanTest, message=gettextRcmdr("The group and block variables must be different."))
        return()
      }
      strLineBlock = paste(" | ", strBlock, " ", sep="") #
      strBlockBlock = paste(", block = ", .activeDataSet, "$",strBlock, sep="")
    } else {
      strLineBlock = "" #
      strBlockBlock = ""
    }
    
    maxpts <- as.integer(tclvalue(asymptoticMaxpts))
    abseps <- as.double(tclvalue(asymptoticAbseps))
    releps <- as.double(tclvalue(asymptoticReleps))
    replications <- as.integer(tclvalue(approximateReplications))
    closeDialog()
    strExactOption = ""
    subset <- tclvalue(subsetVariable) #
    subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) "" else paste(", subset=", subset, sep="") #
    test <- as.character(tclvalue(testVariable))
    
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("tapply(", paste(.activeDataSet,"$", response, sep=""),
                       ", ", paste(.activeDataSet,"$", group, sep=""), ", median, na.rm=TRUE)", sep=""))
    
    if (test == "default"){
      doItAndPrint(paste("friedman_test(", response, " ~ ", group, strLineBlock, ", data=", .activeDataSet, subset, ")", sep=""))           
    } else if (test == "asympt"){
      doItAndPrint(paste("friedman_test(", response, " ~ ", group, strLineBlock, ", data=", .activeDataSet, subset, 
                         ", distribution=asymptotic(maxpts = ", maxpts, ", abseps = ", abseps, ", releps = ", releps, " ) )", sep=""))
    } else {
      doItAndPrint(paste("friedman_test(", response, " ~ ", group, strLineBlock, ", data=", .activeDataSet, subset, 
                         ", distribution=approximate(B = ", replications, " ) )", sep=""))
    }
    
    #doItAndPrint(paste("friedman_test(", response, " ~ ", group, block, ", data=", .activeDataSet, subset, ")", sep=""))
    
    # code from coin package examples:
    ### Wilcoxon-Nemenyi-McDonald-Thompson test
    ### Hollander & Wolfe (1999), page 295
    #if (require("multcomp")) {
      ### all pairwise comparisons
      justDoIt(paste("objMultTest <- symmetry_test(",response ," ~ ",group, strLineBlock, ", data = ",.activeDataSet, ",
			teststat = 'max',
                     xtrafo = function(data)
                     trafo(data, factor_trafo = function(x)
                     model.matrix(~ x - 1) %*% t(contrMat(table(x), 'Tukey'))
                     ),
                     ytrafo = function(data)
                     trafo(data, numeric_trafo = rank", strBlockBlock, ")
      )", sep=""))
			### a global test, again
      justDoIt("objGlobalPValue = pvalue(objMultTest)")
      logger(paste("Global p-value:", objGlobalPValue, sep=""))
      ### simultaneous P-values for all pair comparisons
      logger("Pairwise comparisons of groups:")
      doItAndPrint("pvalue(objMultTest, method = 'single-step')")
    #}
    
    
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="friedman_test")
  
  optionsFrame <- tkframe(top) #
  
  radioButtons(optionsFrame, name="test", buttons=c("default", "asympt", "approx"), 
               labels=c("Default", "Asymptotic", "Approximate"), 
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
  
  
  pairwiseVariable <- tclVar("0")
  pairwiseCheckBox <- tkcheckbutton(optionsFrame, variable=pairwiseVariable)
  
  subsetBox() #
  
  tkgrid(getFrame(groupBox), labelRcmdr(variablesFrame, text="    "), getFrame(responseBox), getFrame(blockBox), sticky="nw") # 
  tkgrid(variablesFrame, sticky="nw") #
  
  tkgrid(testFrame, sticky="nw") #
  tkgrid(optionsFrame, sticky="w")
  
  tkgrid(labelRcmdr(distributionsFrame, text="Distribution options", fg="blue"),sticky="w") #
  
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
  
  
  tkgrid(labelRcmdr(optionsFrame, text="Pairwise comparisons of groups"), pairwiseCheckBox, sticky="w")
  
  
  tkgrid(subsetFrame, sticky="w") #
  tkgrid(buttonsFrame, sticky="w") #
  dialogSuffix(rows=6, columns=2) #
}  
