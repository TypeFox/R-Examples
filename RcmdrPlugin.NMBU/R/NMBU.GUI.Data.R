# $Id: NMBU.GUI.Data.R 35 2014-01-10 21:17:26Z khliland $

##
## GUI functions for the Data menu
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

########################
# Mean center variables
meanCenter <- function(){
  initializeDialog(title=gettextRcmdr("Mean center variables"))
  xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variables (pick one or more)"),
                          selectmode="multiple")
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    if (length(x) == 0) {
      errorCondition(recall=meanCenter, message=gettextRcmdr("You must select one or more variables."))
      return()
    }
    xx <- paste('"', x, '"', sep="")
    .activeDataSet <- ActiveDataSet()
    command <- paste("scale(", .activeDataSet, "[,c(", paste(xx, collapse=","),
                     ")], scale=FALSE)", sep="")
    result <- justDoIt(command)
#    assign(".Z", result, envir=.GlobalEnv)
#    logger(paste(".Z <- ", command, sep=""))
    doItAndPrint(paste(".Z <- ", command, sep=""))
    for (i in 1:length(x)){
      Z <- paste("Z.", x[i], sep="")
      if (is.element(Z, Variables())) {
        if ("no" == tclvalue(checkReplace(Z))){
          if (GrabFocus()) tkgrab.release(top)
          tkdestroy(top)
          next
        }
      }
      justDoIt(paste(.activeDataSet, "$", Z, " <- .Z[,", i, "]", sep=""))
      logger(paste(.activeDataSet, "$", Z, " <- .Z[,", i, "]", sep=""))
    }
    remove(.Z, envir=.GlobalEnv)
    logger("remove(.Z)")
    if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="scale")
  tkgrid(getFrame(xBox), sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=2, columns=1)
}

#####################################
# Sort data based on one variable
sortData <- function(){
  initializeDialog(title=gettextRcmdr("Sort data"))
  variablesFrame <- tkframe(top)
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Variable to use for sorting (pick one)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    if (length(x) != 1){
      errorCondition(recall=sortData, message=gettextRcmdr("You must select one variable."))
      return()
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    justDoIt(paste(.activeDataSet, " <- ", .activeDataSet, "[order(", .activeDataSet, "$", x, "),]", sep=""))
    logger(paste(.activeDataSet, " <- ", .activeDataSet, "[order(", .activeDataSet, "$", x, "),]", sep=""))
    
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="order")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="nw")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=2, columns=1)
}

###################################
## Create data sequence
createSequence <- function(){
  initializeDialog(title=gettextRcmdr("Create data sequence"))
  
  # Boks for signifikans ved ny variabel
  dataSetName   <- tclVar("Dataset")
  dataSetFrame  <- tkframe(top)
  variableName  <- tclVar("Sequence")
  variableFrame <- tkframe(top)
  fromName  <- tclVar("0")
  fromFrame <- tkframe(top)
  toName  <- tclVar("1")
  toFrame <- tkframe(top)
  byName  <- tclVar("0.1")
  byFrame <- ttklabelframe(top,text=gettextRcmdr("Spacing/sequence length"))
  repsName  <- tclVar("1")
  repsFrame <- ttklabelframe(top,text=gettextRcmdr("Repeats"))
  fullFrame <- tkframe(top)
  dataSet   <- ttkentry(dataSetFrame, width="15", textvariable=dataSetName)
  variable  <- ttkentry(variableFrame, width="15", textvariable=variableName)
  from <- ttkentry(fromFrame, width="6", textvariable=fromName)
  to   <- ttkentry(toFrame, width="6", textvariable=toName)
  by   <- ttkentry(byFrame, width="6", textvariable=byName)
  reps <- ttkentry(repsFrame, width="6", textvariable=repsName)
  onOK <- function(){ # Actions to perform
    the.dataSet <- tclvalue(dataSetName)
    the.variable <- tclvalue(variableName)
    the.from <- tclvalue(fromName)
    the.to <- tclvalue(toName)
    the.by <- tclvalue(byName)
    the.reps <- tclvalue(repsName)
    type <- as.character(tclvalue(byOutVariable))
    typeReps <- as.character(tclvalue(repsOutVariable))
    closeDialog()
    if(trim.blanks(the.dataSet) == gettextRcmdr("") || trim.blanks(the.from) == gettextRcmdr("") || trim.blanks(the.to) == gettextRcmdr("") || trim.blanks(the.by) == gettextRcmdr("")){
      errorCondition(recall=createSequence, message=gettextRcmdr("Please fill in all text fields"))
      return()
    }
	if(exists(the.dataSet)){
		command <- paste(the.dataSet, "$", the.variable, " <- rep(seq(", the.from, ", ", the.to, ", ", ifelse(type=="by", "by", "length.out"), "=", the.by, "),",ifelse(typeReps=="each","each=",""),the.reps,")", sep="")
	} else {
		command <- paste(the.dataSet, " <- data.frame(", the.variable, " = rep(seq(", the.from, ", ", the.to, ", ", ifelse(type=="by", "by", "length.out"), "=", the.by, "),",ifelse(typeReps=="each","each=",""),the.reps,"))", sep="")
	}
    doItAndPrint(command)
    activeDataSet(the.dataSet)
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="seq", model=TRUE)
  tkgrid(labelRcmdr(dataSetFrame, text=gettextRcmdr("Name of data set:")), dataSet, sticky="w")
  tkgrid(dataSetFrame, sticky="w", row=1, column=1, columnspan=2)
  tkgrid(labelRcmdr(variableFrame, text=gettextRcmdr("Name of variable:")), variable, sticky="w")
  tkgrid(variableFrame, sticky="w", row=2, column=1, columnspan=2)
  tkgrid(labelRcmdr(fromFrame, text=gettextRcmdr("From:")), from, sticky="w")
  tkgrid(fromFrame, sticky="w", row=3, column=1)
  tkgrid(labelRcmdr(toFrame, text=gettextRcmdr("To:    ")), to, sticky="w")
  tkgrid(toFrame, sticky="nw", row=4, column=1)
  tkgrid(labelRcmdr(byFrame, text=gettextRcmdr("Value")), by, sticky="w")
  radioButtonsNMBU(byFrame,name="byOut", buttons=c("by", "length.out"), values=c("by", "length.out"), initialValue = "by",
                  labels=gettextRcmdr(c("By (spacing)", "Length.out (sequence length)")))#, title=gettextRcmdr("Type"))
  tkgrid(byFrame, sticky="w", row=3, column=2,rowspan=2)
  tkgrid(labelRcmdr(repsFrame, text=gettextRcmdr("Times")), reps, sticky="w")
  radioButtonsNMBU(repsFrame,name="repsOut", buttons=c("all", "each"), values=c("all", "each"), initialValue = "all",
                  labels=gettextRcmdr(c("All (sequence)", "Each (elements)")))#, title=gettextRcmdr("Type"))
  tkgrid(repsFrame, sticky="w", row=5, column=2,rowspan=2)
  tkgrid(buttonsFrame, stick="s", row=7, column=1, columnspan=2)
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=7, columns=2)
}


#################################
# New data set for OS X users
newDataSetOSX <- function() {
  initializeDialog(title=gettextRcmdr("New Data Set"))
  dsname <- tclVar(gettextRcmdr("Dataset"))
  entryDsname <- ttkentry(top, width="20", textvariable=dsname)
  onOK <- function(){
    dsnameValue <- trim.blanks(tclvalue(dsname))
    if (dsnameValue == "") {
      errorCondition(recall=newDataSetOSX,
                     message=gettextRcmdr("You must enter the name of a data set."))
      return()
    }
    if (!is.valid.name(dsnameValue)) {
      errorCondition(recall=newDataSetOSX,
                     message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
      return()
    }
    if (is.element(dsnameValue, listDataSets())) {
      if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
        newDataSetOSX()
        return()
      }
    }
    command <- "edit(data.frame(var1=0))"
    result <- justDoIt(command)
    result <- as.data.frame(lapply(result, function(x) if (is.character(x)) factor(x) else x))
    if (class(result)[1] !=  "try-error"){ 
#			assign(dsnameValue, result, envir=.GlobalEnv)
# 			logger(paste(dsnameValue, "<-", command))
	  doItAndPrint(paste(dsnameValue, "<-", command))
      if (nrow(get(dsnameValue)) == 0){
        #        	if (eval(parse(text=paste("nrow(", dsnameValue, ")"))) == 0){
        errorCondition(recall=newDataSet, message=gettextRcmdr("empty data set."))
        return()
        }
      activeDataSet(dsnameValue)
      }
    closeDialog()
    tkfocus(CommanderWindow())
    }
  OKCancelHelp(helpSubject="edit.data.frame")
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  tkgrid(buttonsFrame, columnspan="2", sticky="w")
  tkgrid.configure(entryDsname, sticky="w")
  dialogSuffix(rows=2, columns=2, focus=entryDsname)
  }

#################################################
# Automatic import of data
auto.import.GUI <- function(){
  the.import <- NULL
  w <- options("warn")$warn
  options(warn=-1)
  try(the.import <- auto.import(), silent=TRUE)
  options(warn=w)
  if(is.null(the.import)){
    Message(message=gettextRcmdr("Automatic import failed"),
            type="error")
    return()
  }
  initializeDialog(title=gettextRcmdr("Automatic import from clipboard"))
  dsname <- tclVar(gettextRcmdr("Dataset"))
  entryDsname <- ttkentry(top, width="20", textvariable=dsname)
  onOK <- function(){
    dsnameValue <- trim.blanks(tclvalue(dsname))
    if(is.null(the.import)){
      errorCondition(recall=auto.import.GUI, message=gettextRcmdr("Automatic import failed."))
    }
    if (dsnameValue == "") {
      errorCondition(recall=auto.import.GUI,
                     message=gettextRcmdr("You must enter a name for the data set."))
      return()
    }
    if (!is.valid.name(dsnameValue)) {
      errorCondition(recall=auto.import.GUI,
                     message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
      return()
    }
    if (is.element(dsnameValue, listDataSets())) {
      if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
        auto.import.GUI()
        return()
      }
    }
    MAC <- Sys.info()[1] == "Darwin"
    cl <- ifelse(MAC,"pipe(\'pbpaste\')", "\'clipboard\'")
    command <- paste(dsnameValue, " <- read.table(", cl, ", strip.white=",the.import$strip.white,", sep='",the.import$sep,"', na.strings='NA', header=",the.import$header,", dec='",the.import$dec,"')", sep="")
#    assign(dsnameValue, the.import[[1]], envir=.GlobalEnv)
#    logger(command)
    doItAndPrint(command)
    if (nrow(get(dsnameValue)) == 0){
      #        	if (eval(parse(text=paste("nrow(", dsnameValue, ")"))) == 0){
      errorCondition(recall=auto.import.GUI, message=gettextRcmdr("empty data set."))
      return()
      } else {
        doItAndPrint(paste("str(",dsnameValue,")",sep=""))
      }
    activeDataSet(dsnameValue)
    closeDialog()
    tkfocus(CommanderWindow())
    }
  OKCancelHelp(helpSubject="edit.data.frame")
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
  tkgrid(buttonsFrame, columnspan="2", sticky="w")
  tkgrid.configure(entryDsname, sticky="w")
  dialogSuffix(rows=2, columns=2, focus=entryDsname)
  }


#################################################
# Delete active data set
deleteActiveDataSet <- function(){
  eval(parse(text=paste('rm(',"'",ActiveDataSet(),"',pos=1)",sep="")))
  ActiveDataSet(ActiveDataSet())
  
  models <- listAllModels()
  for(i in 1:length(models)){
    dataSet <- as.character(get(models[i])$call$data)
    if (!exists(dataSet,where=1)){
      eval(parse(text=paste('rm(',"'",models[i],"',pos=1)",sep="")))    
    }
  }
}


########################
# Update factor
updateFactor <- function(){
  initializeDialog(title=gettextRcmdr("Update factor (remove orphaned levels)"))
  xBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factors (pick one or more)"),
                          selectmode="multiple")
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    if (length(x) == 0) {
      errorCondition(recall=updateFactor, message=gettextRcmdr("You must select one or more factors."))
      return()
    }
    xx <- paste('"', x, '"', sep="")
    .activeDataSet <- ActiveDataSet()
	for(i in 1:length(x)){
		command <- paste(.activeDataSet, "[,'", x[i],"'] <- factor(", .activeDataSet, "[,'", x[i],"'])", sep="")
		justDoIt(command)
		logger(command)
	}
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="scale")
  tkgrid(getFrame(xBox), sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=2, columns=1)
}


########################
# Export last output to clipboard in Excel friendly format
onExport  <- function(){
	doItAndPrint("export.instructions()")}
export.instructions <- function(){
	cat(paste("Clipboard export in Excel friendly format is achieved by pressing",
				"<Ctrl+e> (decimal point .) or",
				"<Ctrl+E> (decimal point ,) immediately after producing the wanted",
				"model or model output. To re-run a code line in the Script window",
				"place the cursor on it and press <Ctrl+r>.\n", sep="\n"))
	nothing <- popOutput()
}
onExportE <- function() doItAndPrint("toClipboard(popOutput(), dec='.', rcmdr=TRUE)")
onExportN <- function() doItAndPrint("toClipboard(popOutput(), dec=',', rcmdr=TRUE)")
toClipboard <- function(object, dec=".", sep="\t", quote=FALSE, rcmdr=FALSE, ...){
	x <- try(xtable(object),silent = TRUE)
	if(class(x)[1] == "xtable"){
		df <- as.data.frame(x)
		write.table(df,file="clipboard", sep=sep, dec=dec, quote=quote, ...)
		cat("Successfully copied data to the clipboard.\nCheck alignment of headers after pasting.\n")
	} else {
		if(rcmdr){
			stop("Try re-running the line producing wanted model/model output (Ctrl+r) before export (Ctrl+e/Ctrl+E).")
		} else {
			stop("Unsuccessful export. Supported classes among: methods(xtable)")
		}
	}
}

