#.onAttach <- function(libname, pkgname){
#    resetlfoptions() ####Ev. alte Einstellungen laden!
#    if (!interactive()) return()
#    Rcmdr <- options()$Rcmdr
#    plugins <- Rcmdr$plugins
#    if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
#        Rcmdr$plugins <- c(plugins, pkgname)
#        options(Rcmdr=Rcmdr)
#        closeCommander(ask=FALSE, ask.save=TRUE)
#        Commander()
#        }
#    }

.onAttach <- function(libname, pkgname){
    resetlfoptions()
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

if(getRversion() >= "2.15.1"){
     utils::globalVariables(c("top",
                              "buttonsFrame",
                              "Recode",
                              "readDataSet",
                              "periodVariable",
                              "periodFrame",
                              "methodVariable",
                              "methodFrame",
                              "styleVariable",
                              "styleFrame",
                              "thresbreaksVariable",
                              "thresbreaksFrame",
                              "linearRegressionModel",
                              "unitVariable",
                              "unitFrame",
                              "poolingVariable",
                              "tableVariable",
                              "dsnameValue",
                              "poolingFrame",
                              "tableFrame"), add = TRUE)}

readlfdatasheet <- function(){
initializeDialog(title=gettextRcmdr("New Low Flow Data"))
optionsFrame <- tkframe(top)
dsname <- tclVar(gettextRcmdr("LFdata"))
entrylfname <- ttkentry(optionsFrame, width="20", textvariable=dsname)
radioButtons(optionsFrame,
             "style",
             buttons=c("GRDC", "HZB","LFU","TU"), 
             labels=gettextRcmdr(c("GRDC", "HZB","LfU-Bayern","TU")),
             title=gettextRcmdr("Type of data sheet"),
             initialValue=getOption("RcmdrPlugin.lfstat")$datasheet)
baseflowVariable <- tclVar(getOption("RcmdrPlugin.lfstat")$baseflow)
baseflowCheckBox <- tkcheckbutton(optionsFrame, variable=baseflowVariable)
hyinit <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$hyear)))
entryhyear <- ttkentry(optionsFrame, width = "2",textvariable=hyinit)

 onOK <- function(){
		closeDialog()

                file <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}}')))
                style <-tclvalue(styleVariable)
                options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(datasheet = style)))

                bf <- tclvalue(baseflowVariable) == "1"
                options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(baseflow = bf)))
                hydro <-as.numeric(tclvalue(hyinit))
                options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hyear = hydro)))                
                dsnameValue <- trim.blanks(tclvalue(dsname))
                if (dsnameValue == ""){
			errorCondition(recall=readDataSet,
				message=gettextRcmdr("You must enter a name for the data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=readDataSet,
				message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
			return()
		}
                if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
				readDataSet()
				return()
			}
		}

command <- paste('readlfdata("', file,'", type="',style,'",hyearstart =',hydro,',baseflow =',bf,')',sep = "")

#Changes as J.Fox in 2013 version
doItAndPrint(paste(dsnameValue, " <- ", command, sep=""))
activeDataSet(dsnameValue)
                
#result <- justDoIt(command)
#if (class(result)[1] !=  "try-error"){
#     gassign(dsnameValue, result)
#     activeDataSet(dsnameValue)
#   }
tkfocus(CommanderWindow())
}
OKCancelHelp(helpSubject="readlfdata")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter name for low flow data:")), entrylfname, sticky="w")
tkgrid(optionsFrame, sticky="w")
tkgrid(styleFrame, sticky="w")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Starting month of the hydrological year:")), entryhyear, sticky="w")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Calculate Baseflow")), baseflowCheckBox, sticky="w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=5, columns=2)
}

#########################
#Create LF Data         #
#########################

createlfdatacalc <- function(){
initializeDialog(title = gettextRcmdr("Convert data to lfobj"))
dateFrame <- tkframe(top)
optionsFrame <- tkframe(top)
nameFrame <-  tkframe(top)
dsname <- tclVar(gettextRcmdr("LFdata"))
entrylfname <- ttkentry(nameFrame, width="20", textvariable=dsname)

variablesBox <- variableListBox(top, Variables(), title=gettextRcmdr("Variable containing flow data:"))

  startdate <- tclVar(gettextRcmdr("dd/mm/yyyy"))
  entrydate <- ttkentry(dateFrame, width="12", textvariable=startdate)
  baseflowVariable <- tclVar(getOption("RcmdrPlugin.lfstat")$baseflow)
  baseflowCheckBox <- tkcheckbutton(optionsFrame, variable=baseflowVariable)
  hyinit <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$hyear)))
  entryhyear <- ttkentry(optionsFrame, width = "2",textvariable=hyinit)

onOK <- function(){
  flowVar <- getSelection(variablesBox)
  closeDialog()
  start <- tclvalue(startdate)
  hyear <- tclvalue(hyinit)
  bf <- tclvalue(baseflowVariable) == "1"
  options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(baseflow = bf)))
              
  options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hyear = hyear)))                

  if(length(flowVar) == 0){
    errorCondition(recall=Recode, message=gettextRcmdr("You must select a variable."))
			return()
		}
    dsnameValue <- trim.blanks(tclvalue(dsname))
                if (dsnameValue == ""){
			errorCondition(recall=readDataSet,
				message=gettextRcmdr("You must enter a name for the data set."))
			return()
		}
		if (!is.valid.name(dsnameValue)){
			errorCondition(recall=readDataSet,
				message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
			return()
		}
                if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
				readDataSet()
				return()
			}
		}

  
  command <- paste("createlfobj(x =ts(",ActiveDataSet(),"$",flowVar,") , startdate = \"",start,"\",hyearstart =",hyear,",baseflow =",bf,")",sep = "")
  doItAndPrint(paste(dsnameValue, " <- ", command, sep=""))
#  result <- justDoIt(command)

#if (class(result)[1] !=  "try-error"){
#     gassign(dsnameValue, result)
#     activeDataSet(dsnameValue)
#   }
tkfocus(CommanderWindow())
  
}#End onOK
OKCancelHelp(helpSubject="createlfobj")

tkgrid(labelRcmdr(nameFrame, text=gettextRcmdr("Enter name for low flow data:")), entrylfname, sticky="w")
tkgrid(nameFrame, sticky="w")
tkgrid(getFrame(variablesBox), sticky="nw")
tkgrid(labelRcmdr(dateFrame, text = gettextRcmdr("Start date of the series")),entrydate,
               sticky = "w")
tkgrid(dateFrame, sticky = "w")

tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Starting month of the hydrological year:")), entryhyear, sticky="w")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Calculate Baseflow")), baseflowCheckBox, sticky="w")
tkgrid(optionsFrame, sticky = "w")

tkgrid(buttonsFrame,sticky = "w")
dialogSuffix(rows=6, columns=2)
}

#########################
# Update LF data        #
#########################

updatelfcalc  <- function(){
initializeDialog(title=gettextRcmdr("Update Flow Data"))
optionsFrame <- tkframe(top)
baseflowVariable <- tclVar(getOption("RcmdrPlugin.lfstat")$baseflow)
baseflowCheckBox <- tkcheckbutton(optionsFrame, variable=baseflowVariable)
hyinit <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$hyear)))
entryhyear <- ttkentry(optionsFrame, width = "2",textvariable=hyinit)

 onOK <- function(){
		closeDialog()

                bf <- tclvalue(baseflowVariable) == "1"
                options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(baseflow = bf)))
                hydro <-as.numeric(tclvalue(hyinit))
                options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hyear = hydro)))                
command <- paste('assign("',ActiveDataSet(),'",createlfobj(', ActiveDataSet(),', hyearstart =',hydro,',baseflow =',bf,'))',sep = "")
result <- doItAndPrint(command)
activeDataSet(ActiveDataSet())
tkfocus(CommanderWindow())
}
OKCancelHelp(helpSubject="createlfobj")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Starting month of the hydrological year:")), entryhyear, sticky="w")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Calculate Baseflow")), baseflowCheckBox, sticky="w")
tkgrid(optionsFrame, sticky = "w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=5, columns=2)
}
