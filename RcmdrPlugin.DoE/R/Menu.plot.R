Menu.plot <- function(){
    ## view and change response names
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))

    putRcmdr("curplotfac", character(0))
    putRcmdr("potentialplotfac", names(factor.names(eval(parse(text=.activeDataSet)))))
    
   onSelect <- function(){
     if (tclvalue(tcl(other.factors$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- potentialplotfac[as.numeric(unlist(strsplit(tclvalue(tcl(other.factors$listbox, "curselection")), " ")))+1]
     putRcmdr("potentialplotfac", setdiff(potentialplotfac,add))
     putRcmdr("curplotfac", c(getRcmdr("curplotfac"),add))
     add <- NULL
     tcl(other.factors$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialplotfac,collapse=" ")))
     other.factors$varlist <- potentialplotfac
     tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curplotfac,collapse=" ")))
     sel.factors$varlist <- curplotfac
   }
   onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(sel.factors$listbox, "curselection"))=="") return()
     add <- curplotfac[as.numeric(unlist(strsplit(tclvalue(tcl(sel.factors$listbox, "curselection")), " ")))+1]
     putRcmdr("curplotfac", setdiff(curplotfac,add))
     putRcmdr("potentialplotfac", c(getRcmdr("potentialplotfac"),add))
     add <- NULL
     tcl(sel.factors$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialplotfac,collapse=" ")))
     other.factors$varlist <- potentialplotfac
     tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curplotfac,collapse=" ")))
     sel.factors$varlist <- curplotfac
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
      hilf <- curplotfac
      if (!is.null(response.names(get(ActiveDataSet())))){
          command <- paste("temp <- ", ActiveDataSet(), "; response.names(temp) <- NULL;")
          if (length(hilf)==0)
             command <- paste(command, "plot(temp); rm(temp)")
          else
             command <- paste(command, "plot(temp, select = c(", 
                   paste(dquote(hilf),collapse=","),")); rm(temp)")
      }
      else{   
      if (length(hilf)==0) 
         command <- paste("plot(",ActiveDataSet(),")")
      else 
         command <- paste("plot(",ActiveDataSet(),", select = c(", paste(dquote(hilf),collapse=","),"))")
      }
      doItAndPrint(command)
      closeDialog(top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
 }

   initializeDialog(title=gettextRcmdr("Plot design for selected factors ..."))
    selFrame <- ttklabelframe(top, text=gettextRcmdr("Select or unselect factors for plotting"))
    estbuttonFrame <- ttkframe(selFrame)
    selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"),
            foreground = "darkgreen", command = onSelect,
            default = "normal", borderwidth = 3)
    tkgrid(selectButton)
    deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"),
            foreground = "darkgreen", command = onDeselect,
            default = "normal", borderwidth = 3)
    tkgrid(deselectButton)

    putRcmdr("sel.factors", variableListBox(selFrame, variableList=curplotfac, listHeight=10, 
        title="Current plot factors",selectmode="multiple"))
    putRcmdr("other.factors", variableListBox(selFrame, variableList=potentialplotfac, listHeight=10, 
        title="Potential further plot factors",selectmode="multiple"))

         tkconfigure(sel.factors$listbox, listvariable=tclVar(paste(curplotfac,collapse=" ")))
         sel.factors$varlist <- curplotfac
         tkconfigure(other.factors$listbox, listvariable=tclVar(paste(potentialplotfac,collapse=" ")))
         other.factors$varlist <- potentialplotfac
         tkgrid(getFrame(other.factors), estbuttonFrame, getFrame(sel.factors))
         tkgrid(selFrame,sticky="n")

    OKCancelHelp(helpSubject="Menu.plot")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=3)
}