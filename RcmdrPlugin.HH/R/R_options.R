R_options <- function() {

  initializeDialog(title=gettextRcmdr("R Options"))

  parFrame <- tkframe(top) 
  pchVar <- tclVar(gettextRcmdr("16"))
  pchEntry <- tkentry(parFrame, width=3, textvariable=pchVar)      

  onOK <- function(){

    pch <- tclvalue(pchVar)
    if ("" == pch) pch="16"
    pch <- trim.blanks(pch)

    old.warn <- options(warn=-1)
    pchValue <- as.numeric(pch)
    if (is.na(pchValue)) {
      if (nchar(pch) > 1) pchValue <- substring(pch,1,1)
      pchValue <- paste('"', pchValue, '"', sep="")
    }
    else
      pchValue <- pch
    options(old.warn)

    command <- paste('trellis.par.set(list(
                 superpose.symbol=list(pch=rep(', pchValue, ', length(trellis.par.get("superpose.symbol")$pch))),
                 plot.symbol=list(pch=', pchValue, ')))', sep='')
    doItAndPrint(command)
    
    command <- paste('par(pch=',pchValue,')', sep='')
    doItAndPrint(command)

 
    stars <- (tclvalue(starsVariable) ==  "1")
    command <- paste('options(show.signif.stars=', stars, ')')
    doItAndPrint(command)
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="options")
 
  checkBoxes(frame="starsFrame", boxes=c("stars"), initialValues=c("1"),
             labels=gettextRcmdr(c("Show significance stars")))
  
  tkgrid(tklabel(parFrame, text=gettextRcmdr("Plotting characters")), pchEntry, stick="w")
  tkgrid(parFrame, sticky="w")
  tkgrid(starsFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix(rows=2, columns=1)
}

## source("~/HH-R.package/RcmdrPlugin.HH/R/R_options.R")
