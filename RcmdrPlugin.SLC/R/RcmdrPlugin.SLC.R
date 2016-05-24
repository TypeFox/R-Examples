# Some Rcmdr dialogs for the SLC package

.onAttach <- function(libname, pkgname){
  if (!interactive()) return()
  Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins
  if (!pkgname %in% plugins) {
    Rcmdr$plugins <- c(plugins, pkgname)
    options(Rcmdr=Rcmdr)
    if("package:Rcmdr" %in% search()) {
      if(!getRcmdr("autoRestart")) {
        closeCommander(ask=FALSE, ask.save=TRUE)
        Commander()
      }
    }
    else {
      Commander()
    }
  }
}

if (getRversion() >= '2.15.1') globalVariables(c('top','namefile','datafile','tkbutton',
                                                 'buttonsFrame','length.baseline'))

Rcmdr.slcestimates <- function(){
require(SLC)

getfile <- function() {
  command <- "tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} 
                {{Data files} {.dat}} {{All files} *}'))"
  gassign("namefile",justDoIt(command))
  if (namefile == "") return();
  result <- justDoIt("data.frame(values=scan(namefile))")
  if (class(result)[1] != "try-error"){
    gassign("datafile", result)
    activeDataSet("datafile")
  }
}

initializeDialog(title=gettextRcmdr("Slope and Level Change Estimates"))
lengthVar <- tclVar("1")
lengthEntry <- tkentry(top, width="12", textvariable=lengthVar)

onOK <- function(){
  closeDialog()

        if (is.numeric(datafile[,1]) == FALSE){
           errorCondition(recall=Rcmdr.slcestimates, message="Invalid Data Type: Original data must be numeric.")
            return()
            }

        length.baseline <- as.numeric(tclvalue(lengthVar))
#        command <- paste("as.numeric(",tclvalue(lengthVar),")", sep="")
#        gassign("length.baseline", justDoIt(command))
        if ( (is.numeric(length.baseline) == FALSE) | (length.baseline < 1) ) {
            errorCondition(recall=Rcmdr.slcestimates, message="The length of baseline phase must be greater than 0.")
            return()
            }
	tkfocus(CommanderWindow())

        doItAndPrint(paste("slcestimates(datafile$values,",length.baseline,")",sep=''))      
	tkfocus(CommanderWindow())
}

OKCancelHelp(helpSubject="slcestimates")
tkgrid(labelRcmdr(top, text=""))
tkgrid(labelRcmdr(top, text=gettextRcmdr(" - Original Data should be loaded before OK - "), fg="blue"), sticky="w")
tkgrid(labelRcmdr(top, text=""), sticky="w")
tkgrid(tkbutton(top, text="Select Data File", command=getfile))
tkgrid(labelRcmdr(top, text=""), sticky="w")
tkgrid(labelRcmdr(top, text=gettextRcmdr(" - Baseline Phase Length Specification - "), fg="blue"),sticky="w")
tkgrid(labelRcmdr(top, text=""), sticky="w")
tkgrid(tklabel(top, text="Baseline Phase Length:"), lengthEntry)
tkgrid(labelRcmdr(top, text=""), sticky="w")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
tkgrid.configure(lengthEntry, sticky="w")
dialogSuffix(rows=10, columns=1)
}

Rcmdr.help.SLC <- function(){
   require(SLC)
   doItAndPrint("help(\"SLC\")")
   invisible(NULL)
}

Rcmdr.help.RcmdrPlugin.SLC <- function(){
   require(SLC)
   doItAndPrint("help(\"RcmdrPlugin.SLC\")")
   invisible(NULL)
}

