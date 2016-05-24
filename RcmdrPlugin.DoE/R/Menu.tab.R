Menu.tab <- function(){
    ## view and change response names
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))

    putRcmdr("curtabfac", character(0))
    putRcmdr("potentialtabfac", names(factor.names(eval(parse(text=.activeDataSet)))))
    
   onSelect <- function(){
     if (tclvalue(tcl(other.factors$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- potentialtabfac[as.numeric(unlist(strsplit(tclvalue(tcl(other.factors$listbox, "curselection")), " ")))+1]
     putRcmdr("potentialtabfac", setdiff(potentialtabfac,add))
     putRcmdr("curtabfac", c(getRcmdr("curtabfac"),add))
     add <- NULL
     tcl(other.factors$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialtabfac,collapse=" ")))
     other.factors$varlist <- potentialtabfac
     tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curtabfac,collapse=" ")))
     sel.factors$varlist <- curtabfac
   }
   onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(sel.factors$listbox, "curselection"))=="") return()
     add <- curtabfac[as.numeric(unlist(strsplit(tclvalue(tcl(sel.factors$listbox, "curselection")), " ")))+1]
     putRcmdr("curtabfac", setdiff(curtabfac,add))
     putRcmdr("potentialtabfac", c(getRcmdr("potentialtabfac"),add))
     add <- NULL
     tcl(sel.factors$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialtabfac,collapse=" ")))
     other.factors$varlist <- potentialtabfac
     tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curtabfac,collapse=" ")))
     sel.factors$varlist <- curtabfac
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
      hilf <- curtabfac
      if (length(hilf)==0) 
         command <- paste("table(",ActiveDataSet(),")")
      else 
         command <- paste("table(",ActiveDataSet(),"[, c(", paste(dquote(hilf),collapse=","),")])")
      doItAndPrint(command)
      closeDialog(top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
 }

   initializeDialog(title=gettextRcmdr("Tabulate design for selected factors ..."))
    selFrame <- ttklabelframe(top, text=gettextRcmdr("Select or unselect factors for tabulating"))
    estbuttonFrame <- ttkframe(selFrame)
    selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"),
            foreground = "darkgreen", command = onSelect,
            default = "normal", borderwidth = 3)
    tkgrid(selectButton)
    deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"),
            foreground = "darkgreen", command = onDeselect,
            default = "normal", borderwidth = 3)
    tkgrid(deselectButton)

    putRcmdr("sel.factors", variableListBox(selFrame, variableList=curtabfac, listHeight=10, 
        title="Current table factors",selectmode="multiple"))
    putRcmdr("other.factors", variableListBox(selFrame, variableList=potentialtabfac, listHeight=10, 
        title="Potential further table factors",selectmode="multiple"))

         tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curtabfac,collapse=" ")))
         sel.factors$varlist <- curtabfac
         tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialtabfac,collapse=" ")))
         other.factors$varlist <- potentialtabfac
         tkgrid(getFrame(other.factors), estbuttonFrame, getFrame(sel.factors))
         tkgrid(selFrame,sticky="n")

    OKCancelHelp(helpSubject="Menu.tab")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=3)
}