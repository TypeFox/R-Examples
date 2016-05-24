Menu.responses <- function(){
    ## view and change response names
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))

    putRcmdr("curresp", response.names(eval(parse(text=.activeDataSet))))
    putRcmdr("potentialresp", setdiff(listNumeric(),c(curresp,names(di$factor.names))))

   onSelect <- function(){
     if (tclvalue(tcl(other.numerics$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- potentialresp[as.numeric(unlist(strsplit(tclvalue(tcl(other.numerics$listbox, "curselection")), " ")))+1]
     putRcmdr("potentialresp", setdiff(potentialresp,add))
     putRcmdr("curresp", c(getRcmdr("curresp"),add))
     add <- NULL
     tcl(other.numerics$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.numerics$listbox, listvariable=tclVar(paste(potentialresp,collapse=" ")))
     other.numerics$varlist <- potentialresp
     tkconfigure(sel.responses$listbox, listvariable=tclVar(paste(curresp,collapse=" ")))
     sel.responses$varlist <- curresp
   }
   onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(sel.responses$listbox, "curselection"))=="") return()
     add <- curresp[as.numeric(unlist(strsplit(tclvalue(tcl(sel.responses$listbox, "curselection")), " ")))+1]
     putRcmdr("curresp", setdiff(curresp,add))
     putRcmdr("potentialresp", c(getRcmdr("potentialresp"),add))
     add <- NULL
     tcl(sel.responses$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.numerics$listbox, listvariable=tclVar(paste(potentialresp,collapse=" ")))
     other.numerics$varlist <- potentialresp
     tkconfigure(sel.responses$listbox, listvariable=tclVar(paste(curresp,collapse=" ")))
     sel.responses$varlist <- curresp
 }
dquote <- function(obj){
    ## quote vector elements for use as character vector in a command
    aus <- rep("",length(obj))
    wopt <- options("warn")[[1]]
    options(warn=-1)
    for (i in 1:length(obj)) if (is.na(as.numeric(obj[i]))) {
            if (length(grep('"',obj[i])>0))
            aus[i] <- paste("'",obj[i],"'",sep="") 
            else
            aus[i] <- paste('"',obj[i],'"',sep="") 
            }
          else aus[i] <- obj[i]
    options(warn=wopt)
    aus
}

onOK <- function(){
      hilf <- curresp
      if (length(hilf)==0) 
      command <- paste("response.names(",ActiveDataSet(),") <- NULL")
      else
      command <- paste("response.names(",ActiveDataSet(),") <- c(", paste(dquote(hilf),collapse=","),")")
      doItAndPrint(command)
      putRcmdr(".activeDataSet", ActiveDataSet())
      activeDataSet(.activeDataSet)
      closeDialog(top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
 }

   initializeDialog(title=gettextRcmdr("Modify response settings for the design ..."))
    selFrame <- ttklabelframe(top, text=gettextRcmdr("Select or unselect responses"))
    estbuttonFrame <- ttkframe(selFrame)
    selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"),
            foreground = "darkgreen", command = onSelect,
            default = "normal", borderwidth = 3)
    tkgrid(selectButton)
    deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"),
            foreground = "darkgreen", command = onDeselect,
            default = "normal", borderwidth = 3)
    tkgrid(deselectButton)

    putRcmdr("sel.responses", variableListBox(selFrame, variableList=curresp, listHeight=10, 
        title="Current responses",selectmode="multiple"))
    putRcmdr("other.numerics", variableListBox(selFrame, variableList=potentialresp, listHeight=10, 
        title="Potential further responses",selectmode="multiple"))

         tkconfigure(sel.responses$listbox, listvariable=tclVar(paste(curresp,collapse=" ")))
         sel.responses$varlist <- curresp
         tkconfigure(other.numerics$listbox, listvariable=tclVar(paste(potentialresp,collapse=" ")))
         other.numerics$varlist <- potentialresp
         tkgrid(getFrame(other.numerics), estbuttonFrame, getFrame(sel.responses))
         tkgrid(selFrame,sticky="n")

    OKCancelHelp(helpSubject="Menu.responses")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=3)
}