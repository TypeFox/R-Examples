multitablecalc <- function(){
initializeDialog(title = gettextRcmdr("Calculate Indices of various stations"))
 optionsFrame <- tkframe(top)
radioButtons(optionsFrame,
             "method",
            buttons=c("MRC", "IRS"), 
             labels=gettextRcmdr(c("MRC", "IRS")),
             title=gettextRcmdr("Recession Method"),
              initialValue = getOption("RcmdrPlugin.lfstat")$recessionmethod)
lfobjBox1 <- variableListBox(top, listlfobj(),
		selectmode="multiple", title=gettextRcmdr("Select low flow objects: "))
indices <- c("meanflow", "Q95", "MAM1", "MAM7", "MAM10", "MAM30", "MAM90", "baseflowindex", "recession")
indicesBox <- variableListBox(top,
                            indices,
                            title=gettextRcmdr("Indices"),
                              selectmode="multiple",
                            initialSelection=0:(length(indices)-1))

  seglen <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$segLength))
  entryseglength <- ttkentry(optionsFrame, width="2", textvariable=seglen)
  thrlevel <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$threslevelrec))
  entrythrlevel <- ttkentry(optionsFrame, width = "3", textvariable = thrlevel)
  tmean <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$tmean))
  entrytmean <- ttkentry(optionsFrame, width= "3", textvariable = tmean)
onOK <-  function(){
  lfobj1 <- getSelection(lfobjBox1)
  ind <- getSelection(indicesBox)
  closeDialog()
method <- tclvalue(methodVariable)
seg <- tclvalue(seglen)
thr <- tclvalue(thrlevel)
trim <- tclvalue(tmean)
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(threslevelrec = thr)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(recessionmethod = method)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(segLength = seg)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(tmean = trim)))
  indname2 <- NULL
  lfname <- NULL
  
  for(ii in seq_along(lfobj1)){
    lfname <- paste(lfname,lfobj1[ii],sep = if(ii == 1){""}else{","})
  }
  for (ii in seq_along(ind)){
    indname2 <- paste(indname2,ind[ii],sep =if(ii == 1){"\""}else{"\",\""})
  }
  
  command <- paste('multistationsreport(',lfname,',indices =c(',indname2,'"),recessionmethod = "',method,'",recessionseglength = ',seg,',recessionthreshold = ',thr,', recessiontrimIRS = ',trim,')',sep="")
  doItAndPrint(command)
  tkfocus(CommanderWindow())
  } #End on OK


OKCancelHelp(helpSubject="multistationsreport")
	tkgrid(getFrame(lfobjBox1),getFrame(indicesBox), sticky="nw")

tkgrid(optionsFrame,sticky="w")
tkgrid(methodFrame,sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Recession threshold level:")), entrythrlevel, sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Recession segment length:")), entryseglength, sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Recession trim level (IRS only):")), entrytmean, sticky = "w")
tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows=2, columns=2)
}
