ChisqHypothesesPlot <- function(){
    initializeDialog(title=gettextRcmdr("Chisq Distribution"))
    degfreeVar <- tclVar("1")
    degfreeEntry <- tkentry(top, width="6", textvariable=degfreeVar)
    alphaLeftVar <- tclVar("")
    alphaLeftEntry <- tkentry(top, width="6", textvariable=alphaLeftVar)
    alphaRightVar <- tclVar(".05")
    alphaRightEntry <- tkentry(top, width="6", textvariable=alphaRightVar)
    obsValueVar <- tclVar("")
    obsValueEntry <- tkentry(top, width="6", textvariable=obsValueVar)
    xmaxVar <- tclVar("")
    xmaxEntry <- tkentry(top, width="6", textvariable=xmaxVar)
    ymaxVar <- tclVar("")
    ymaxEntry <- tkentry(top, width="6", textvariable=ymaxVar)
    onOK <- function(){
      closeDialog()
      degfree <- as.numeric(tclvalue(degfreeVar))
      alphaLeft <- as.numeric(tclvalue(alphaLeftVar))
      alphaRight <- as.numeric(tclvalue(alphaRightVar))
      obsValue <- as.numeric(tclvalue(obsValueVar))
      xmax <- as.numeric(tclvalue(xmaxVar))
      ymax <- as.numeric(tclvalue(ymaxVar))
      
      command <- "old.par <- par(oma=c(1,1,0,5))"
      justDoIt(command)
      logger(command)
      
      command.xlim <- if (!is.na(xmax)) 
        paste(", xlim=c(0, ", xmax, ")", sep="")
      else ""
      
      command.ylim <- if (!is.na(ymax)) 
        paste(", ylim=c(0, ", ymax, ")", sep="")
      else ""
      
      command <- "chisq.setup("
      if (!is.na(degfree))
        command <- paste(command, "df=", degfree, sep="")
      else
        degfree <-1

      command <- paste(command,
                       command.xlim,
                       command.ylim,
                       sep="")

      command <- paste(command, ")", sep="")
      justDoIt(command)
      logger(command)

      command.alpha <- ", alpha=c("
      if (!is.na(alphaLeft)) command.alpha <-
        paste(command.alpha, alphaLeft, ", ", sep="")
      command.alpha <-  paste(command.alpha, alphaRight, ")", sep="")

      shade <- if (!is.na(alphaLeft)) "outside" else "right"

      
      command <- paste("chisq.curve(",
                       "df=", degfree,
                       command.alpha, 
                       ", shade='", shade,
                       "', col='blue')", sep="")
      doItAndPrint(command)

      ## Observed Value
      if (!is.na(obsValue)) {
        command <- paste("chisq.observed(", obsValue, ", col='green'",
                       ", df=", degfree,
                       ")",  sep="")
        doItAndPrint(command)
      }

      command <- "par(old.par)"
      justDoIt(command)
      logger(command)
      
      tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="dnorm")
    tkgrid(tklabel(top, text=gettextRcmdr("df (degrees of freedom)")), degfreeEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("left alpha")), alphaLeftEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("right alpha")), alphaRightEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("Observed Value")), obsValueEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("Chisq max (right-hand side)")), xmaxEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("f(Chisq) max (top side)")), ymaxEntry, sticky="e")
    tkgrid(buttonsFrame, columnspan=1, sticky="w")
    tkgrid.configure(degfreeEntry, sticky="w")
    tkgrid.configure(alphaLeftEntry, sticky="w")
    tkgrid.configure(alphaRightEntry, sticky="w")
    tkgrid.configure(obsValueEntry, sticky="w")
    tkgrid.configure(xmaxEntry, sticky="w")
    tkgrid.configure(ymaxEntry, sticky="w")
    dialogSuffix(rows=8, columns=1, focus=degfreeEntry)
    }

## source("~/HH-R.package/RcmdrPlugin.HH/R/ChisqHypothesesPlot.R")
## source("~/HH-R.package/HH/R/chisq.curve.R")
