.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                options(Rcmdr=Rcmdr)
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
            Commander()
        }
    }
}

########################################################################################################
##########                                Visual analysis                                     ##########
########################################################################################################

############################
# Graphical representation #
############################

Rcmdr.graph <- function(){

  initializeDialog(title=c("Graphical display"))
  
  designBox <- variableListBox(top, variableList=c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  onOK <- function(){

    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.graph, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    dat <- tclvalue(dataFileVariable)
    if (dat == "active") {
      .activeDataSet <- ActiveDataSet()
      if (length(ActiveDataSet()) == 0){
        errorCondition(recall=Rcmdr.graph, message=c("You must select an active data set."))
        return()
      }
      command <- paste("graph(design = ",design, ", data = ", .activeDataSet, ")", sep="")
    }
    if (dat == "load") {
    
      if (exists("dataFile") == 0){
        errorCondition(recall=Rcmdr.graph, message=c("You must select a data file."))
        return()
      }
      command <- paste("graph(design = ",design, ", data = read.table(dataFile))", sep="")
    }

    closeDialog()    
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="graph")
  
  datFrame <- tkframe(top)

  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))

  buttonFrame <- tkframe(datFrame)

  tkgrid(getFrame(designBox), sticky="we")
  tkgrid(labelRcmdr(top, text=""))
 
  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))

  tkgrid(buttonsFrame)
  dialogSuffix()
}

####################################
# Plot measure of central tendency #
####################################

Rcmdr.graph.CL <- function(){

  initializeDialog(title=c("Display central location"))

  designBox <- variableListBox(top, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))

  clFrame <- tkframe(top)

  clBox <- variableListBox(clFrame, variableList= c("Mean", "Median", "Broadened median", "Trimmed mean", "M-estimator"), listHeight=4, title=c("Select the measure of central tendency"))

  onOK <- function(){

    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.graph.CL, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}

    cl <- getSelection(clBox)
    if (length(cl) == 0){
      errorCondition(recall=Rcmdr.graph.CL, message=c("You must select a measure of central tendency."))
      return()
    }
    if (cl == "Mean") {CL <- paste("\"mean\"")}
    if (cl == "Median") {CL <- paste("\"median\"")}
    if (cl == "Broadened median") {CL <- paste("\"bmed\"")}
    if (cl == "Trimmed mean") {CL <- paste("\"trimmean\"")}
    if (cl == "M-estimator") {CL <- paste("\"mest\"")}

    trim <- tclvalue(trimVariable)
    if (cl == "Trimmed Mean"){
    if (trim == ""){
        errorCondition(recall=Rcmdr.graph.CL, message=c("You must indicate the % of observations to be removed from the end of the distribution."))
        return()
      }
    }
    trim <- as.numeric(trim)

    mest <- tclvalue(mestVariable)
    if (cl == "M-estimator"){
    if (mest == ""){
        errorCondition(recall=Rcmdr.graph.CL, message=c("You must indicate the value for the constant K."))
        return()
      }
    }
    mest <- as.numeric(mest)
        
    dat <- tclvalue(dataFileVariable)

    if (cl == "Mean" | cl == "Median" | cl == "Broadened median"){
      if (dat == "active") {
          .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", data = ", .activeDataSet, ")", sep="")
      }
      if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", data = read.table(dataFile))", sep="")
      }
    }
    
    if (cl == "Trimmed mean"){
      if (dat == "active") {
          .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", tr = ",trim, ", data = ", .activeDataSet, ")", sep="")
      }
      if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", tr = ",trim, ", data = read.table(dataFile))", sep="")
      }
    }

    if (cl == "M-estimator"){
      if (dat == "active") {
          .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", tr = ",mest, ", data = ", .activeDataSet, ")", sep="")
      }
      if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.CL, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.CL(design = ",design, ", CL = ",CL, ", tr = ",mest, ", data = read.table(dataFile))", sep="")
      }
    }

    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="graph.CL")
  
  constantFrame<- tkframe(clFrame)

  trimFrame <- tkframe(constantFrame)
  trimVariable <- tclVar("0.2")
  trimField <- ttkentry(trimFrame, width="6", textvariable=trimVariable)

  mestFrame <- tkframe(constantFrame)
  mestVariable <- tclVar("1.28")
  mestField <- ttkentry(mestFrame, width="6", textvariable=mestVariable)
  
  datFrame <- tkframe(top)

  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))

  buttonFrame <- tkframe(datFrame)

  tkgrid(getFrame(designBox), sticky="nw")
  tkgrid(labelRcmdr(top, text=""))

  tkgrid(labelRcmdr(trimFrame, text=c("For trimmed mean:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(trimFrame, text=c("% of observations to be removed: ")), trimField, sticky="w")
  tkgrid(trimField, sticky="w")
  tkgrid(labelRcmdr(constantFrame, text=""))
  tkgrid(labelRcmdr(mestFrame, text=c("For M-estimator:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(mestFrame, text=c("Value for the constant K: ")), mestField, sticky="w")
  tkgrid(mestField, sticky="w")
  tkgrid(trimFrame, sticky = "nw")
  tkgrid(mestFrame, sticky = "nw")
  tkgrid(getFrame(clBox), labelRcmdr(clFrame, text= "      "), constantFrame, sticky = "nw")
  tkgrid(clFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))

  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))

  tkgrid(buttonsFrame)
  dialogSuffix()
}

################################
# Plot estimate of variability #
################################

Rcmdr.graph.VAR <- function(){
  
  initializeDialog(title=c("Display variability"))

  selectionFrame <- tkframe(top)

  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))

  varBox <- variableListBox(selectionFrame, variableList= c("Range lines", "Range bars", "Trended range"), listHeight=4, title=c("Select the measure of variability."))

  clFrame <- tkframe(top)

  clBox <- variableListBox(clFrame, variableList= c("Mean", "Median", "Broadened median", "Trimmed mean", "M-estimator"), listHeight=4, title=c("Select the measure of central tendency"))

  onOK <- function(){

    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}

    var <- getSelection(varBox)
    if (length(var) == 0){
      errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a measure of variability."))
      return()
    }
    if (var == "Range lines") {VAR <- paste("\"RL\"")}
    if (var == "Range bars") {VAR <- paste("\"RB\"")}
    if (var == "Trended range") {VAR <- paste("\"TR\"")}

    cl <- getSelection(clBox)
    if (length(cl) == 0){
      errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a measure of central tendency."))
      return()
    }
    if (cl == "Mean") {CL <- paste("\"mean\"")}
    if (cl == "Median") {CL <- paste("\"median\"")}
    if (cl == "Broadened median") {CL <- paste("\"bmed\"")}
    if (cl == "Trimmed mean") {CL <- paste("\"trimmean\"")}
    if (cl == "M-estimator") {CL <- paste("\"mest\"")}

    trim <- tclvalue(trimVariable)
    if (cl == "Trimmed Mean"){
      if (trim == ""){
        errorCondition(recall=Rcmdr.graph.VAR, message=c("You must indicate the % of observations to be removed from the end of the distribution."))
        return()
      }
    }
    trim <- as.numeric(trim)

    mest <- tclvalue(mestVariable)
    if (cl == "M-estimator"){
      if (mest == ""){
        errorCondition(recall=Rcmdr.graph.VAR, message=c("You must indicate the value for the constant K."))
        return()
      }
    }
    mest <- as.numeric(mest)

    dat <- tclvalue(dataFileVariable)

    dataSet <- tclvalue(dataSetVariable)
    if (dataSet == 0){dataset <- paste("\"regular\"")}
    if (dataSet == 1){dataset <- paste("\"trimmed\"")}
       
    if (var == "Range lines" | var == "Trended range"){
      if (dat == "active") {
	.activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select an active data set."))
          return()
        }
        command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", data = ", .activeDataSet, ")", sep="")
      }
      if (dat == "load") {
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a data file."))
          return()
        }
        command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", data = read.table(dataFile))", sep="")
      }
    } 

    if (var == "Range bars"){
      if (cl == "Mean" | cl == "Median" | cl == "Broadened median"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", data = read.table(dataFile))", sep="")
        }
      }
      if (cl == "Trimmed mean"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", tr = ",trim, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", tr = ",trim, ", data = read.table(dataFile))", sep="")
        }
      }
      if (cl == "M-estimator"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", tr = ",mest, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.VAR, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.VAR(design = ",design, ", VAR = ",VAR, ", dataset = ",dataset, ", CL = ",CL, ", tr = ",mest, ", data = read.table(dataFile))", sep="")
        }
      }
    }   
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="graph.VAR")

  constantFrame<- tkframe(clFrame)

  trimFrame <- tkframe(constantFrame)
  trimVariable <- tclVar("0.2")
  trimField <- ttkentry(trimFrame, width="6", textvariable=trimVariable)

  mestFrame <- tkframe(constantFrame)
  mestVariable <- tclVar("1.28")
  mestField <- ttkentry(mestFrame, width="6", textvariable=mestVariable)
  
  datFrame <- tkframe(top)

  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))

  buttonFrame <- tkframe(datFrame)

  optionsFrame <- tkframe(top)
  dataSetVariable <- tclVar("0")
  dataSetCheckBox <- tkcheckbutton(optionsFrame, variable=dataSetVariable)

  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(varBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
 
  tkgrid(labelRcmdr(top, text=c("                                          For range bars:"), fg="blue"), sticky = "w")

  tkgrid(labelRcmdr(trimFrame, text=c("For trimmed mean:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(trimFrame, text=c("% of observations to be removed: ")), trimField, sticky="w")
  tkgrid(trimField, sticky="w")
  tkgrid(labelRcmdr(constantFrame, text=""))
  tkgrid(labelRcmdr(mestFrame, text=c("For M-estimator:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(mestFrame, text=c("Value for the constant K: ")), mestField, sticky="w")
  tkgrid(mestField, sticky="w")
  tkgrid(trimFrame, sticky = "nw")
  tkgrid(mestFrame, sticky = "nw")
  tkgrid(getFrame(clBox), labelRcmdr(clFrame, text= "      "), constantFrame, sticky = "nw")
  tkgrid(clFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
 
  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(optionsFrame, text=c("Remove the extreme values")), dataSetCheckBox, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)
  dialogSuffix()
}

##########################
# Plot estimate of trend #
##########################

Rcmdr.graph.TREND <- function(){
  
  initializeDialog(title=c("Display trend"))

  selectionFrame <- tkframe(top)

  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))

  trendBox <- variableListBox(selectionFrame, variableList= c("Vertical line plot", "Trend lines (Least Squares regression)", "Trend lines (Split-middle)", "Trend lines (Resistant trend line fitting)", "Running medians (batch size 3)", "Running medians (batch size 5)", "Running medians (batch size 4 averaged by pairs)"), listHeight=4, title=c("Select the trend visualization"))

  clFrame <- tkframe(top)

  clBox <- variableListBox(clFrame, variableList= c("Mean", "Median", "Broadened median", "Trimmed mean", "M-estimator"), listHeight=4, title=c("Select the measure of central tendency"))

  onOK <- function(){

    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}

    trend <- getSelection(trendBox)
    if (length(trend) == 0){
      errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a trend visualization."))
      return()
    }
    if (trend == "Vertical line plot") {TREND <- paste("\"VLP\"")}
    if (trend == "Trend lines (Least Squares regression)") {TREND <- paste("\"LSR\"")}
    if (trend == "Trend lines (Split-middle)") {TREND <- paste("\"SM\"")}
    if (trend == "Trend lines (Resistant trend line fitting)") {TREND <- paste("\"RTL\"")}
    if (trend == "Running medians (batch size 3)") {TREND <- paste("\"RM3\"")}
    if (trend == "Running medians (batch size 5)") {TREND <- paste("\"RM5\"")}
    if (trend == "Running medians (batch size 4 averaged by pairs)") {TREND <- paste("\"RM42\"")}

    cl <- getSelection(clBox)
    if (length(cl) == 0){
      errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a measure of central tendency."))
      return()
    }
    if (cl == "Mean") {CL <- paste("\"mean\"")}
    if (cl == "Median") {CL <- paste("\"median\"")}
    if (cl == "Broadened median") {CL <- paste("\"bmed\"")}
    if (cl == "Trimmed mean") {CL <- paste("\"trimmean\"")}
    if (cl == "M-estimator") {CL <- paste("\"mest\"")}

    trim <- tclvalue(trimVariable)
    if (cl == "Trimmed Mean"){
      if (trim == ""){
        errorCondition(recall=Rcmdr.graph.TREND, message=c("You must indicate the % of observations to be removed from the end of the distribution."))
        return()
      }
    }
    trim <- as.numeric(trim)

    mest <- tclvalue(mestVariable)
    if (cl == "M-estimator"){
      if (mest == ""){
        errorCondition(recall=Rcmdr.graph.TREND, message=c("You must indicate the value for the constant K."))
        return()
      }
    }
    mest <- as.numeric(mest)

    dat <- tclvalue(dataFileVariable)

    if (trend == "Trend lines (Least Squares regression)" | trend == "Trend lines (Split-middle)" | trend == "Trend lines (Resistant trend line fitting)" | trend == "Running medians (batch size 3)" | trend == "Running medians (batch size 5)" | trend == "Running medians (batch size 4 averaged by pairs)"){
      if (dat == "active") {
	.activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select an active data set."))
          return()
        }
        command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", data = ", .activeDataSet, ")", sep="")
      }
      if (dat == "load") {
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a data file."))
          return()
        }
        command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", data = read.table(dataFile))", sep="")	
      }
    }

    if (trend == "Vertical line plot"){
      if (cl == "Mean" | cl == "Median" | cl == "Broadened median"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", data = read.table(dataFile))", sep="")	
        }

      }
      if (cl == "Trimmed mean"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", tr = ",trim, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", tr = ",trim, ", data = read.table(dataFile))", sep="")	
        }
      }
      if (cl == "M-estimator"){
        if (dat == "active") {
	  .activeDataSet <- ActiveDataSet()
          if (length(ActiveDataSet()) == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select an active data set."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", tr = ",mest, ", data = ", .activeDataSet, ")", sep="")
        }
        if (dat == "load") {
          if (exists("dataFile") == 0){
            errorCondition(recall=Rcmdr.graph.TREND, message=c("You must select a data file."))
            return()
          }
          command <- paste("graph.TREND(design = ",design, ", TREND = ",TREND, ", CL = ",CL, ", tr = ",mest, ", data = read.table(dataFile))", sep="")	
        }
      }
    }
      
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="graph.TREND")

  constantFrame<- tkframe(clFrame)

  trimFrame <- tkframe(constantFrame)
  trimVariable <- tclVar("0.2")
  trimField <- ttkentry(trimFrame, width="6", textvariable=trimVariable)

  mestFrame <- tkframe(constantFrame)
  mestVariable <- tclVar("1.28")
  mestField <- ttkentry(mestFrame, width="6", textvariable=mestVariable)
  
  datFrame <- tkframe(top)

  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))

  buttonFrame <- tkframe(datFrame)

  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(trendBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
 
  tkgrid(labelRcmdr(top, text=c("                                          For vertical line plot:"), fg="blue"), sticky = "w")

  tkgrid(labelRcmdr(trimFrame, text=c("For trimmed mean:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(trimFrame, text=c("% of observations to be removed: ")), trimField, sticky="w")
  tkgrid(trimField, sticky="w")
  tkgrid(labelRcmdr(constantFrame, text=""))
  tkgrid(labelRcmdr(mestFrame, text=c("For M-estimator:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(mestFrame, text=c("Value for the constant K: ")), mestField, sticky="w")
  tkgrid(mestField, sticky="w")
  tkgrid(trimFrame, sticky = "nw")
  tkgrid(mestFrame, sticky = "nw")
  tkgrid(getFrame(clBox), labelRcmdr(clFrame, text= "      "), constantFrame, sticky = "nw")
  tkgrid(clFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
 
  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))

  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))

  tkgrid(buttonsFrame)
  dialogSuffix()
}

########################################################################################################
##########                            Randomization tests                                     ##########
########################################################################################################

###########################################################
##        Designing single-case experiments              ##
###########################################################

##################################
# Number of possible assignments #
##################################

Rcmdr.quantity <- function(){
  
  initializeDialog(title=c("Number of possible assignments"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.quantity, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    MT <- tclvalue(MTVariable)
    if (des !=  "Multiple Baseline Design") {
      if (MT == ""){
        errorCondition(recall=Rcmdr.quantity, message=c("You must indicate the number of observations."))
        return()
      }
    }
    MT <- as.numeric(MT)
    
    phase <- tclvalue(phaseVariable)
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      if (phase == ""){
        errorCondition(recall=Rcmdr.quantity, message=c("You must indicate the minimum number of observations per phase."))
        return()
      }
    }
    phase <- as.numeric(phase)
    
    atd <- tclvalue(atdVariable)
    if (des == "Alternating Treatments Design"){
      if (atd == ""){
        errorCondition(recall=Rcmdr.quantity, message=c("You must indicate the maximum number of consecutive administrations of the same condition."))
        return()
      }
    }
    atd <- as.numeric(atd)
    
    if (des ==  "Multiple Baseline Design") {
      
      ###
      command <- "tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} {{All files} *}'))"
      starts <- putRcmdr("starts", justDoIt(command)) 
      if (starts == "") return();
      ###
      
      if (exists("starts") == 0){
        errorCondition(recall=Rcmdr.quantity, message=c("You must select a start points file."))
        return()
      }
      starts <- paste("\"", starts, "\"", sep = "")
      command <- paste("quantity(design = ", design, ", starts = ", starts, ")", sep="")
    }
    
    if (des == "Completely Randomized Design" | des == "Randomized Block Design") {
      command <- paste("quantity(design = ", design, ", MT = ", MT, ")", sep="")}
    
    if (des == "Alternating Treatments Design") {command <- paste("quantity(design = ", design, ", MT = ", MT, ", limit = ", atd, ")", sep="")}
    
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design") {command <- paste("quantity(design = ", design, ", MT = ", MT, ", limit = ", phase, ")", sep="")}
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="quantity")
  
  MTFrame <- tkframe(selectionFrame)
  MTVariable <- tclVar("")
  MTField <- ttkentry(MTFrame, width="6", textvariable=MTVariable)
  
  limFrame <- tkframe(top)
  phaseFrame <- tkframe(limFrame)
  phaseVariable <- tclVar("")
  phaseField <- ttkentry(phaseFrame, width="6", textvariable=phaseVariable)
  atdFrame <- tkframe(limFrame)
  atdVariable <- tclVar("")
  atdField <- ttkentry(atdFrame, width="6", textvariable=atdVariable)
  
  startsFrame <- tkframe(top)
  
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=c("Number of observations \n(not necessary for MBD): ")), MTField, sticky="w")
  tkgrid(MTField, sticky="w")
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), MTFrame, sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(phaseFrame, text=c("For phase designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("For alternating treatments designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(phaseFrame, text=c("Minimum number of observations per \nphase: ")), sticky="w")
  tkgrid(phaseField, sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("Maximum number of consecutive administrations \nof the same condition: ")), sticky="w")
  tkgrid(atdField, sticky="w")
  tkgrid(phaseFrame, labelRcmdr(limFrame, text = "         "), atdFrame, sticky = "nw")
  tkgrid(limFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(top, text=c("If the chosen design is a multiple baseline design, a file select window
will appear (after clicking 'OK') where you must select the file containing the possible startpoints."), fg="blue"), sticky = "w")
  tkgrid(startsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)
  dialogSuffix()
  
}

####################################
# Display all possible assignments #
####################################

Rcmdr.assignments <- function(){
  
  initializeDialog(title=c("Display all possible data arrangements"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.assignments, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    MT <- tclvalue(MTVariable)
    if (des !=  "Multiple Baseline Design") {
      if (MT == ""){
        errorCondition(recall=Rcmdr.assignments, message=c("You must indicate the number of observations."))
        return()
      }
    }
    MT <- as.numeric(MT)
    
    phase <- tclvalue(phaseVariable)
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      if (phase == ""){
        errorCondition(recall=Rcmdr.assignments, message=c("You must indicate the minimum number of observations per phase."))
        return()
      }
    }
    phase <- as.numeric(phase)
    
    atd <- tclvalue(atdVariable)
    if (des == "Alternating Treatments Design"){
      if (atd == ""){
        errorCondition(recall=Rcmdr.assignments, message=c("You must indicate the maximum number of consecutive administrations of the same condition."))
        return()
      }
    }
    atd <- as.numeric(atd)
    
    sav <- tclvalue(saveAssignmentsVariable)
    if (sav == "Check") {save <- paste("\"yes\"")}
    
    if (sav == "No") {save <- paste("\"no\"")}
    
    if (des == "Multiple Baseline Design") {
      
      ###
      command <- "tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} {{All files} *}'))"
      starts <- putRcmdr("starts", justDoIt(command))
      if (starts == "") return();
      ###
      
      if (exists("starts") == 0){
        errorCondition(recall=Rcmdr.assignments, message=c("You must select a start points file."))
        return()
      }
      starts <- paste("\"", starts, "\"", sep = "")  
      command <- paste("assignments(design = ", design, ", save = ", save, ", starts = ", starts, ")", sep="")
    }
    
    if (des == "Completely Randomized Design" | des == "Randomized Block Design") {
      command <- paste("assignments(design = ", design, ", MT = ", MT, ", save = ", save, ")", sep="")
    }
    
    if (des == "Alternating Treatments Design") {
      command <- paste("assignments(design = ", design, ", MT = ", MT, ", limit = ", atd, ", save = ", save, ")", sep="")
    }
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design") {
      command <- paste("assignments(design = ", design, ", MT = ", MT, ", limit = ", phase, ", save = ", save, ")", sep="")
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="assignments")
  
  MTFrame <- tkframe(selectionFrame)
  MTVariable <- tclVar("")
  MTField <- ttkentry(MTFrame, width="6", textvariable=MTVariable)
  
  limFrame <- tkframe(top)
  phaseFrame <- tkframe(limFrame)
  phaseVariable <- tclVar("")
  phaseField <- ttkentry(phaseFrame, width="6", textvariable=phaseVariable)
  atdFrame <- tkframe(limFrame)
  atdVariable <- tclVar("")
  atdField <- ttkentry(atdFrame, width="6", textvariable=atdVariable)
  
  startsFrame <- tkframe(top)
  
  savFrame <- tkframe(top)
  saveFrame <- tkframe(savFrame)
  radioButtons(saveFrame, name="saveAssignments", buttons = c("Check", "No"), labels= c(c("YES", "NO")), title=c("Save the possible assignments?"))
  buttonFrame <- tkframe(savFrame)
  
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=c("Number of observations \n(not necessary for MBD): ")), MTField, sticky="w")
  tkgrid(MTField, sticky="w")
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), MTFrame, sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(phaseFrame, text=c("For phase designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("For alternating treatments designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(phaseFrame, text=c("Minimum number of observations per \nphase: ")), sticky="w")
  tkgrid(phaseField, sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("Maximum number of consecutive administrations \nof the same condition: ")), sticky="w")
  tkgrid(atdField, sticky="w")
  tkgrid(phaseFrame, labelRcmdr(limFrame, text = "         "), atdFrame, sticky = "nw")
  tkgrid(limFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(saveAssignmentsFrame, sticky="w")
  tkgrid(labelRcmdr(saveFrame, text=c("Choose 'YES' if you want to save the possible assignments to a text file.
In this case, a file select window will appear (after clicking 'OK').
Create a text file to store the assignments by typing 'filename.txt'
in the window in a folder of your choice. ")), sticky="w")
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(saveFrame, labelRcmdr(savFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(savFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=c("If the chosen design is a multiple baseline design, a file select window will appear (after clicking 'OK')
where you must select the file containing the possible startpoints. If you also chose for the assignments
to be saved, a second window will appear where you must create a text file to store the assignments."), fg="blue"), sticky = "w")
  tkgrid(startsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)
  dialogSuffix()
  
}

################################
# Choose 1 possible assignment #
################################

Rcmdr.selectdesign <- function(){
  
  initializeDialog(title=c("Choose 1 data arrangement"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.selectdesign, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    MT <- tclvalue(MTVariable)
    if (des !=  "Multiple Baseline Design") {
      if (MT == ""){
        errorCondition(recall=Rcmdr.selectdesign, message=c("You must indicate the number of observations."))
        return()
      }
    }
    MT <- as.numeric(MT)
    
    phase <- tclvalue(phaseVariable)
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      if (phase == ""){
        errorCondition(recall=Rcmdr.selectdesign, message=c("You must indicate the minimum number of observations per phase."))
        return()
      }
    }
    phase <- as.numeric(phase)
    
    atd <- tclvalue(atdVariable)
    if (des == "Alternating Treatments Design"){
      if (atd == ""){
        errorCondition(recall=Rcmdr.selectdesign, message=c("You must indicate the maximum number of consecutive administrations of the same condition."))
        return()
      }
    }
    atd <- as.numeric(atd)
    
    if (des == "Multiple Baseline Design"){
      
      command <- paste("selectdesign(design = ", design,")", sep="")
      
    }
    if (des == "Completely Randomized Design" | des == "Randomized Block Design"){
      command <- paste("selectdesign(design = ", design, ", MT = ", MT, ")", sep="")
    }  
    if (des == "Alternating Treatments Design"){
      command <- paste("selectdesign(design = ", design, ", MT = ", MT, ", limit = ", atd, ")", sep="")
    }
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      command <- paste("selectdesign(design = ", design, ", MT = ", MT, ", limit = ", phase, ")", sep="")
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
    
  }
  
  OKCancelHelp(helpSubject="selectdesign")
  
  MTFrame <- tkframe(selectionFrame)
  MTVariable <- tclVar("")
  MTField <- ttkentry(MTFrame, width="6", textvariable=MTVariable)
  
  limFrame <- tkframe(top)
  phaseFrame <- tkframe(limFrame)
  phaseVariable <- tclVar("")
  phaseField <- ttkentry(phaseFrame, width="6", textvariable=phaseVariable)
  atdFrame <- tkframe(limFrame)
  atdVariable <- tclVar("")
  atdField <- ttkentry(atdFrame, width="6", textvariable=atdVariable)
  
  startsFrame <- tkframe(top)
  
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=""))
  tkgrid(labelRcmdr(MTFrame, text=c("Number of observations \n(not necessary for MBD): ")), MTField, sticky="w")
  tkgrid(MTField, sticky="w")
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), MTFrame, sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(phaseFrame, text=c("For phase designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("For alternating treatments designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(phaseFrame, text=c("Minimum number of observations per \nphase: ")), sticky="w")
  tkgrid(phaseField, sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("Maximum number of consecutive administrations \nof the same condition: ")), sticky="w")
  tkgrid(atdField, sticky="w")
  tkgrid(phaseFrame, labelRcmdr(limFrame, text = "         "), atdFrame, sticky = "nw")
  tkgrid(limFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(top, text=c("If the chosen design is a multiple baseline design, a file select window will appear (after clicking 'OK')
where you must select the file containing the possible startpoints."), fg="blue"), sticky = "w")
  tkgrid(startsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)
  dialogSuffix()
  
}

##########################################
##               ANALYSIS               ##
##########################################

###########################
# Observed test statistic #
###########################

Rcmdr.observed <- function(){
  
  initializeDialog(title=c("Observed test statistic"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  statisticBox <- variableListBox(selectionFrame, variableList= c("A-B", "B-A", "|A-B|", "PA-PB (only for ABA and ABAB)", "PB-PA (only for ABA and ABAB)", "|PA-PB| (only for ABA and ABAB)", "AA-BB (only for ABA and ABAB)", "BB-AA (only for ABA and ABAB)", "|AA-BB| (only for ABA and ABAB)"), listHeight=4, title=c("Select the test statistic"))
  
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.observed, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    stat <- getSelection(statisticBox)
    if (length(stat) == 0){
      errorCondition(recall=Rcmdr.observed, message=c("You must select a test statistic."))
      return()
    }
    if (stat == "A-B") {statistic <- paste("\"A-B\"")}
    if (stat == "B-A") {statistic <- paste("\"B-A\"")}
    if (stat == "|A-B|") {statistic <- paste("\"|A-B|\"")}
    if (stat == "PA-PB (only for ABA and ABAB)") {statistic <- paste("\"PA-PB\"")}
    if (stat == "PB-PA (only for ABA and ABAB)") {statistic <- paste("\"PB-PA\"")}
    if (stat == "|PA-PB| (only for ABA and ABAB)") {statistic <- paste("\"|PA-PB|\"")}
    if (stat == "AA-BB (only for ABA and ABAB)") {statistic <- paste("\"AA-BB\"")}
    if (stat == "BB-AA (only for ABA and ABAB)") {statistic <- paste("\"BB-AA\"")}
    if (stat == "|AA-BB| (only for ABA and ABAB)") {statistic <- paste("\"|AA-BB|\"")}
    
    dat <- tclvalue(dataFileVariable)
    if (dat == "active") {
      .activeDataSet <- ActiveDataSet()
      if (length(ActiveDataSet()) == 0){
        errorCondition(recall=Rcmdr.observed, message=c("You must select an active data set."))
        return()
      }
      command <- paste("observed(design = ",design, ", statistic = ",statistic, ", data = ", .activeDataSet, ")", sep="")
    }
    if (dat == "load") {
      if (exists("dataFile") == 0){
        errorCondition(recall=Rcmdr.observed, message=c("You must select a data file."))
        return()
      }
      command <- paste("observed(design = ",design, ", statistic = ",statistic, ", data = read.table(dataFile))", sep="")
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="observed")
  
  datFrame <- tkframe(top)
  
  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))
  
  buttonFrame <- tkframe(datFrame)
  
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(statisticBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(dataFileFrame, sticky="w")
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)  
  dialogSuffix()
}

##############################
# Randomization distribution #
##############################

Rcmdr.distribution <- function(){
  
  initializeDialog(title=c("Randomization distribution"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  statisticBox <- variableListBox(selectionFrame, variableList= c("A-B", "B-A", "|A-B|", "PA-PB (only for ABA and ABAB)", "PB-PA (only for ABA and ABAB)", "|PA-PB| (only for ABA and ABAB)", "AA-BB (only for ABA and ABAB)", "BB-AA (only for ABA and ABAB)", "|AA-BB| (only for ABA and ABAB)"), listHeight=4, title=c("Select the test statistic"))
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.distribution, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    stat <- getSelection(statisticBox)
    if (length(stat) == 0){
      errorCondition(recall=Rcmdr.distribution, message=c("You must select a test statistic."))
      return()
    }
    if (stat == "A-B") {statistic <- paste("\"A-B\"")}
    if (stat == "B-A") {statistic <- paste("\"B-A\"")}
    if (stat == "|A-B|") {statistic <- paste("\"|A-B|\"")}
    if (stat == "PA-PB (only for ABA and ABAB)") {statistic <- paste("\"PA-PB\"")}
    if (stat == "PB-PA (only for ABA and ABAB)") {statistic <- paste("\"PB-PA\"")}
    if (stat == "|PA-PB| (only for ABA and ABAB)") {statistic <- paste("\"|PA-PB|\"")}
    if (stat == "AA-BB (only for ABA and ABAB)") {statistic <- paste("\"AA-BB\"")}
    if (stat == "BB-AA (only for ABA and ABAB)") {statistic <- paste("\"BB-AA\"")}
    if (stat == "|AA-BB| (only for ABA and ABAB)") {statistic <- paste("\"|AA-BB|\"")}
    
    distr <- tclvalue(distributionVariable)
    
    phase <- tclvalue(phaseVariable)
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      if (phase == ""){
        errorCondition(recall=Rcmdr.distribution, message=c("You must indicate the minimum number of observations per phase."))
        return()
      }
    }
    phase <- as.numeric(phase)
    
    atd <- tclvalue(atdVariable)
    if (des == "Alternating Treatments Design"){
      if (atd == ""){
        errorCondition(recall=Rcmdr.distribution, message=c("You must indicate the maximum number of consecutive administrations of the same condition."))
        return()
      }
    }
    atd <- as.numeric(atd)
    
    number <- tclvalue(numberVariable)
    if (distr == "mc"){
      if (number == ""){
        errorCondition(recall=Rcmdr.distribution, message=c("You must indicate the number of randomizations."))
        return()
      }
    }
    number <- as.numeric(number)
    
    sav <- tclvalue(saveAssignmentsVariable)
    if (sav == "Check") {save <- paste("\"yes\"")}
    
    if (sav == "No") {save <- paste("\"no\"")}
    
    dat <- tclvalue(dataFileVariable)
    
    if (des == "Multiple Baseline Design"){
      
      ###
      command <- "tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} {{All files} *}'))"
      starts <- putRcmdr("starts", justDoIt(command))
      if (starts == "") return();
      ###
      
      if (exists("starts") == 0){
        errorCondition(recall=Rcmdr.distribution, message=c("You must select a start points file."))
        return()
      }
      starts <- paste("\"", starts, "\"", sep = "")
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.distribution, message=c("You must select an active data set."))
          return()
        }
        if (distr == "mc"){
          command <- paste("distribution.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", save = ", save, ", data = ", .activeDataSet, ", starts = ", starts, ")", sep="")
        }
        if (distr == "syst"){
          command <- paste("distribution.systematic(design = ", design, ", statistic = ",statistic, ", save = ", save, ", data = ", .activeDataSet, ", starts = ", starts, ")", sep="")
        }
      }
    }
    
    if (des == "Completely Randomized Design" | des == "Randomized Block Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.distribution, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("distribution.systematic(design = ", design, ", statistic = ",statistic, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("distribution.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      
    }
    
    if (des == "Alternating Treatments Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.distribution, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("distribution.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("distribution.random(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", number = ", number, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      
    }
    
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.distribution, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("distribution.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("distribution.random(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", number = ", number, ", save = ", save, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="distribution.systematic")
  
  limFrame <- tkframe(top)
  phaseFrame <- tkframe(limFrame)
  phaseVariable <- tclVar("")
  phaseField <- ttkentry(phaseFrame, width="6", textvariable=phaseVariable)
  atdFrame <- tkframe(limFrame)
  atdVariable <- tclVar("")
  atdField <- ttkentry(atdFrame, width="6", textvariable=atdVariable)
  
  randomFrame <- tkframe(top)
  
  randFrame <- tkframe(randomFrame)
  
  radioButtons(randFrame, name="distribution", buttons= c("syst", "mc"), labels = c(c("Systematic randomization distribution", "Monte Carlo randomization distribution")), title=c("Select the randomization distribution"))
  
  numberFrame <- tkframe(randFrame)
  numberVariable <- tclVar("")
  numberField <- ttkentry(numberFrame, width="6", textvariable=numberVariable)
  
  savFrame <- tkframe(randomFrame)
  saveFrame <- tkframe(savFrame)
  radioButtons(saveFrame, name="saveAssignments", buttons= c("No", "Check"), labels = c(c("NO     ", "YES     ")), title=c("Save the distribution?"))
  savbuttonFrame <- tkframe(savFrame)
  
  datFrame <- tkframe(top)
  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))
  buttonFrame <- tkframe(datFrame)
  
  startsFrame <- tkframe(top)
  
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(statisticBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(phaseFrame, text=c("For phase designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("For alternating treatments designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(phaseFrame, text=c("Minimum number of observations per \nphase: ")), sticky="w")
  tkgrid(phaseField, sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("Maximum number of consecutive administrations \nof the same condition: ")), sticky="w")
  tkgrid(atdField, sticky="w")
  tkgrid(phaseFrame, labelRcmdr(limFrame, text = "         "), atdFrame, sticky = "nw")
  tkgrid(limFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(distributionFrame, sticky = "w")
  tkgrid(labelRcmdr(numberFrame, text=c("    Number of randomizations: ")), numberField, sticky = "w")
  tkgrid(numberFrame, sticky = "w")
  tkgrid(saveAssignmentsFrame, sticky="w")
  tkgrid(labelRcmdr(saveFrame, text=c("Choose 'YES' if you want to save the possible assignments to a text file.
In this case, a file select window will appear (after clicking 'OK').
Create a text file to store the assignments by typing 'filename.txt'
in the window in a folder of your choice.")), sticky="w")
  tkgrid(labelRcmdr(savbuttonFrame, text=""))
  
  tkgrid(saveFrame, labelRcmdr(savFrame, text="  "), savbuttonFrame, sticky="w")
  tkgrid(randFrame, labelRcmdr(randomFrame, text = "      "), savFrame, sticky = "w") 
  tkgrid(randomFrame, sticky = "w") 
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))
  
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(top, text=c("If the chosen design is a multiple baseline design, a file select window will appear (after clicking 'OK')
where you must select the file containing the possible startpoints. If you also chose for the assignments
to be saved, a second window will appear where you must create a text file to store the assignments."), fg="blue"), sticky = "w")
  
  tkgrid(startsFrame, sticky="w")
  
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
  
}

###########
# P value # 
###########

Rcmdr.pvalue <- function(){
  
  initializeDialog(title=c("P-value"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  statisticBox <- variableListBox(selectionFrame, variableList= c("A-B", "B-A", "|A-B|", "PA-PB (only for ABA and ABAB)", "PB-PA (only for ABA and ABAB)", "|PA-PB| (only for ABA and ABAB)", "AA-BB (only for ABA and ABAB)", "BB-AA (only for ABA and ABAB)", "|AA-BB| (only for ABA and ABAB)"), listHeight=4, title=c("Select the test statistic"))
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.pvalue, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    stat <- getSelection(statisticBox)
    if (length(stat) == 0){
      errorCondition(recall=Rcmdr.pvalue, message=c("You must select a test statistic."))
      return()
    }
    if (stat == "A-B") {statistic <- paste("\"A-B\"")}
    if (stat == "B-A") {statistic <- paste("\"B-A\"")}
    if (stat == "|A-B|") {statistic <- paste("\"|A-B|\"")}
    if (stat == "PA-PB (only for ABA and ABAB)") {statistic <- paste("\"PA-PB\"")}
    if (stat == "PB-PA (only for ABA and ABAB)") {statistic <- paste("\"PB-PA\"")}
    if (stat == "|PA-PB| (only for ABA and ABAB)") {statistic <- paste("\"|PA-PB|\"")}
    if (stat == "AA-BB (only for ABA and ABAB)") {statistic <- paste("\"AA-BB\"")}
    if (stat == "BB-AA (only for ABA and ABAB)") {statistic <- paste("\"BB-AA\"")}
    if (stat == "|AA-BB| (only for ABA and ABAB)") {statistic <- paste("\"|AA-BB|\"")}
    
    distr <- tclvalue(distributionVariable)
    
    phase <- tclvalue(phaseVariable)
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design"){
      if (phase == ""){
        errorCondition(recall=Rcmdr.pvalue, message=c("You must indicate the minimum number of observations per phase."))
        return()
      }
    }
    phase <- as.numeric(phase)
    
    atd <- tclvalue(atdVariable)
    if (des == "Alternating Treatments Design"){
      if (atd == ""){
        errorCondition(recall=Rcmdr.pvalue, message=c("You must indicate the maximum number of consecutive administrations of the same condition."))
        return()
      }
    }
    atd <- as.numeric(atd)
    
    number <- tclvalue(numberVariable)
    if (distr == "mc"){
      if (number == ""){
        errorCondition(recall=Rcmdr.pvalue, message=c("You must indicate the number of randomizations."))
        return()
      }
    }
    number <- as.numeric(number)
    
    dat <- tclvalue(dataFileVariable)
    
    if (des == "Multiple Baseline Design"){
      
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select an active data set."))
          return()
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      if (dat == "load") { 
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select a data file."))
          return()
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", data = read.table(dataFile), starts = ", starts, ")", sep="")
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", data = read.table(dataFile), starts = ", starts, ")", sep="")
        }
      }
    }
    
    if (des == "Completely Randomized Design" | des == "Randomized Block Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      if (dat == "load") { 
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select a data file."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", data = read.table(dataFile))", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", number = ", number, ", data = read.table(dataFile))", sep="")
        }
      }
    }
    
    if (des == "Alternating Treatments Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", number = ", number, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      if (dat == "load") { 
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select a data file."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", data = read.table(dataFile))", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", limit = ", atd, ", number = ", number, ",  data = read.table(dataFile))", sep="")
        }
      }
    }
    
    if (des == "AB Phase Design" | des == "ABA Phase Design" | des == "ABAB Phase Design") {
      if (dat == "active") {
        .activeDataSet <- ActiveDataSet()
        if (length(ActiveDataSet()) == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select an active data set."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", data = ", .activeDataSet, ")", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", number = ", number, ", data = ", .activeDataSet, ")", sep="")
        }
      }
      if (dat == "load") { 
        if (exists("dataFile") == 0){
          errorCondition(recall=Rcmdr.pvalue, message=c("You must select a data file."))
          return()
        }
        if (distr == "syst"){
          command <- paste("pvalue.systematic(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", data = read.table(dataFile))", sep="")
        }
        if (distr == "mc"){
          command <- paste("pvalue.random(design = ", design, ", statistic = ",statistic, ", limit = ", phase, ", number = ", number, ",  data = read.table(dataFile))", sep="")
        }
      }
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="pvalue.systematic")
  
  limFrame <- tkframe(top)
  phaseFrame <- tkframe(limFrame)
  phaseVariable <- tclVar("")
  phaseField <- ttkentry(phaseFrame, width="6", textvariable=phaseVariable)
  atdFrame <- tkframe(limFrame)
  atdVariable <- tclVar("")
  atdField <- ttkentry(atdFrame, width="6", textvariable=atdVariable)
  
  randomFrame <- tkframe(top)
  
  radioButtons(randomFrame, name="distribution", buttons= c("syst", "mc"), labels = c(c("Systematic randomization distribution", "Monte Carlo randomization distribution")), title=c("Select the randomization distribution"))
  
  numberFrame <- tkframe(randomFrame)
  numberVariable <- tclVar("")
  numberField <- ttkentry(numberFrame, width="6", textvariable=numberVariable)
  
  datFrame <- tkframe(top)
  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))
  buttonFrame <- tkframe(datFrame)
  
  startsFrame <- tkframe(top)
  
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(statisticBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(phaseFrame, text=c("For phase designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("For alternating treatments designs:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(phaseFrame, text=c("Minimum number of observations per \nphase: ")), sticky="w")
  tkgrid(phaseField, sticky="w")
  tkgrid(labelRcmdr(atdFrame, text=c("Maximum number of consecutive administrations \nof the same condition: ")), sticky="w")
  tkgrid(atdField, sticky="w")
  tkgrid(phaseFrame, labelRcmdr(limFrame, text = "         "), atdFrame, sticky = "nw")
  tkgrid(limFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(numberFrame, text=c("For the Monte Carlo distribution:"), fg="blue"), sticky="w")
  tkgrid(labelRcmdr(numberFrame, text=c("Number of randomizations: ")), numberField, sticky = "w")
  tkgrid(distributionFrame, labelRcmdr(randomFrame, text = "  "), numberFrame, sticky = "w")
  tkgrid(randomFrame, sticky = "w") 
  tkgrid(labelRcmdr(top, text=""))
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(dataFileFrame, sticky="w")
  
  tkgrid(labelRcmdr(buttonFrame, text=""))
  
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(labelRcmdr(top, text=c("If the chosen design is a multiple baseline design, a file select window will appear (after clicking 'OK')
where you must select the file containing the possible startpoints."), fg="blue"), sticky = "w")
  
  tkgrid(startsFrame, sticky="w")
  
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
  
}

########################################################################################################
##########                                  Meta-analysis                                     ##########
########################################################################################################

#############################################
# Effect size measures for single-case data #
#############################################

Rcmdr.ES <- function(){
  
  initializeDialog(title=c("Effect size"))
  
  selectionFrame <- tkframe(top)
  
  designBox <- variableListBox(selectionFrame, variableList= c("AB Phase Design", "ABA Phase Design", "ABAB Phase Design", "Completely Randomized Design", "Alternating Treatments Design", "Randomized Block Design", "Multiple Baseline Design"), listHeight=4, title=c("Select the design type"))
  
  esBox <- variableListBox(selectionFrame, variableList= c("Standardized Mean Difference", "Pooled Standardized Mean Difference", "PND (expected increase)", "PND (expected decrease)", "PEM (expected increase)", "PEM (expected decrease)"), listHeight=4, title=c("Select the effect size measure"))
  
  onOK <- function(){
    
    des <- getSelection(designBox)
    if (length(des) == 0){
      errorCondition(recall=Rcmdr.ES, message=c("You must select a design type."))
      return()
    }
    if (des == "Completely Randomized Design") {design <- paste("\"CRD\"")}
    if (des == "Alternating Treatments Design") {design <- paste("\"ATD\"")}
    if (des == "Randomized Block Design") {design <- paste("\"RBD\"")}
    if (des == "AB Phase Design") {design <- paste("\"AB\"")}
    if (des == "ABA Phase Design") {design <- paste("\"ABA\"")}
    if (des == "ABAB Phase Design") {design <- paste("\"ABAB\"")}
    if (des == "Multiple Baseline Design") {design <- paste("\"MBD\"")}
    
    es <- getSelection(esBox)
    if (length(es) == 0){
      errorCondition(recall=Rcmdr.ES, message=c("You must select an effect size measure."))
      return()
    }
    if (es == "Standardized Mean Difference") {ES <- paste("\"SMD\"")}
    if (es == "Pooled Standardized Mean Difference") {ES <- paste("\"SMDpool\"")}
    if (es == "PND (expected increase)") {ES <- paste("\"PND+\"")}
    if (es == "PND (expected decrease)") {ES <- paste("\"PND-\"")}
    if (es == "PEM (expected increase)") {ES <- paste("\"PEM+\"")}
    if (es == "PEM (expected decrease)") {ES <- paste("\"PEM-\"")}
    
    dat <- tclvalue(dataFileVariable)
    if (dat == "active") {
      .activeDataSet <- ActiveDataSet()
      if (length(ActiveDataSet()) == 0){
        errorCondition(recall=Rcmdr.ES, message=c("You must select an active data set."))
        return()
      }
      command <- paste("ES(design = ",design, ", ES = ",ES, ", data = ", .activeDataSet, ")", sep="")
    }
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="ES")
  
  datFrame <- tkframe(top)
  
  dataFrame <- tkframe(datFrame)
  radioButtons(dataFrame, name="dataFile", buttons = c("active"), labels = c(c("Use the active data set")), title=c("Select the data file"))
  
  buttonFrame <- tkframe(datFrame)
  
  tkgrid(getFrame(designBox), labelRcmdr(selectionFrame, text = "        "), getFrame(esBox), sticky="nw")
  tkgrid(selectionFrame, sticky = "nw")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(dataFileFrame, sticky="w")
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(dataFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)  
  dialogSuffix()
  
}

#########################
# Probability combining #
#########################

Rcmdr.combine <- function(){
  
  initializeDialog(title=c("Combine p-values"))
  
  onOK <- function(){
    
    meth <- tclvalue(combMethodVariable)
    if (meth == "mult") {method <- paste("\"x\"")}
    if (meth == "add") {method <- paste("\"+\"")}
    
    pval <- tclvalue(pvalueFileVariable)
    if (pval == "active") {
      .activeDataSet <- ActiveDataSet()
      if (length(ActiveDataSet()) == 0){
        errorCondition(recall=Rcmdr.combine, message=c("You must select an active data set."))
        return()
      }
      command <- paste("combine(method = ",method, ", pvalues = ", .activeDataSet, ")", sep="")
    }
    
    
    closeDialog()
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  OKCancelHelp(helpSubject="combine")
  
  radioButtons(name = "combMethod", buttons = c("mult", "add"), labels = c(c("Multiplicative", "Additive")), title = c("Select the combining method"))
  
  datFrame <- tkframe(top)
  
  pvalFrame <- tkframe(datFrame)
  radioButtons(pvalFrame, name = "pvalueFile", buttons = c("active"), labels = c(c("Use the active data set")), title = c("Select the p-values file"))
  
  buttonFrame <- tkframe(datFrame)
  
  tkgrid(combMethodFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(pvalueFileFrame, sticky="w")
  tkgrid(labelRcmdr(buttonFrame, text=""))
  tkgrid(pvalFrame, labelRcmdr(datFrame, text="  "), buttonFrame, sticky="w")
  tkgrid(datFrame, sticky = "w")
  tkgrid(labelRcmdr(top, text=""))
  
  tkgrid(buttonsFrame)
  dialogSuffix()
}