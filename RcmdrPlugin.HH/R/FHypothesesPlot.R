FHypothesesPlot <- function(){
    initializeDialog(title=gettextRcmdr("F Distribution"))
    degfree1Var <- tclVar("")
    degfree1Entry <- tkentry(top, width="6", textvariable=degfree1Var)
    degfree2Var <- tclVar("")
    degfree2Entry <- tkentry(top, width="6", textvariable=degfree2Var)
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
      degfree1 <- as.numeric(tclvalue(degfree1Var))
      degfree2 <- as.numeric(tclvalue(degfree2Var))
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
      
      command <- "F.setup("
      if (!is.na(degfree1))
        command <- paste(command, "df1=", degfree1, sep="")
      else
        degfree1 <-1
      if (!is.na(degfree2))
        command <- paste(command, ", df2=", degfree2, sep="")
      else
        degfree2 <- Inf

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

      
      command <- paste("F.curve(",
                       "df1=", degfree1,
                       ", df2=", degfree2,
                       command.alpha, 
                       ", shade='", shade,
                       "', col='blue')", sep="")
      doItAndPrint(command)

      ## Observed Value
      if (!is.na(obsValue)) {
        command <- paste("F.observed(", obsValue, ", col='green'",
                       ", df1=", degfree1,
                       ", df2=", degfree2,
                       ")",  sep="")
        doItAndPrint(command)
      }

      command <- "par(old.par)"
      justDoIt(command)
      logger(command)
      
      tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="dnorm")
    tkgrid(tklabel(top, text=gettextRcmdr("df1 (numerator degrees of freedom)")), degfree1Entry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("df2 (denominator degrees of freedom)")), degfree2Entry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("left alpha")), alphaLeftEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("right alpha")), alphaRightEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("Observed Value")), obsValueEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("F max (right-hand side)")), xmaxEntry, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("f(F) max (top side)")), ymaxEntry, sticky="e")
    tkgrid(buttonsFrame, columnspan=1, sticky="w")
    tkgrid.configure(degfree1Entry, sticky="w")
    tkgrid.configure(degfree2Entry, sticky="w")
    tkgrid.configure(alphaLeftEntry, sticky="w")
    tkgrid.configure(alphaRightEntry, sticky="w")
    tkgrid.configure(obsValueEntry, sticky="w")
    tkgrid.configure(xmaxEntry, sticky="w")
    tkgrid.configure(ymaxEntry, sticky="w")
    dialogSuffix(rows=8, columns=1, focus=degfree1Entry)
    }

## source("~/HH-R.package/RcmdrPlugin.HH/R/FHypothesesPlot.R")
## source("~/HH-R.package/HH/R/F.curve.R")
