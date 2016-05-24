Menu.colremove <- function(){
    ## view and change response names
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))

    putRcmdr("currem", NULL)
    putRcmdr("removable", setdiff(colnames(eval(parse(text=.activeDataSet))),c(names(di$factor.names), di$block.name)))

   onSelect <- function(){
     if (tclvalue(tcl(other.removables$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- removable[as.numeric(unlist(strsplit(tclvalue(tcl(other.removables$listbox, "curselection")), " ")))+1]
     putRcmdr("removable", setdiff(removable,add))
     putRcmdr("currem", c(getRcmdr("currem"),add))
     add <- NULL
     tcl(other.removables$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.removables$listbox, listvariable=tclVar(paste(removable,collapse=" ")))
     other.removables$varlist <- removable
     tkconfigure(sel.columns$listbox, listvariable=tclVar(paste(currem,collapse=" ")))
     sel.columns$varlist <- currem
   }
   onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(sel.columns$listbox, "curselection"))=="") return()
     add <- currem[as.numeric(unlist(strsplit(tclvalue(tcl(sel.columns$listbox, "curselection")), " ")))+1]
     putRcmdr("currem", setdiff(currem,add))
     putRcmdr("removable", c(getRcmdr("removable"),add))
     add <- NULL
     tcl(sel.columns$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.removables$listbox, listvariable=tclVar(paste(removable,collapse=" ")))
     other.removables$varlist <- removable
     tkconfigure(sel.columns$listbox, listvariable=tclVar(paste(currem,collapse=" ")))
     sel.columns$varlist <- currem
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
      hilf <- currem
      if (length(hilf)>0){ 
        command <- paste(ActiveDataSet(), "<- col.remove(",ActiveDataSet(),", c(", paste(dquote(hilf),collapse=","),"))")
        doItAndPrint(command)
        putRcmdr(".activeDataSet", ActiveDataSet())
        activeDataSet(.activeDataSet)
      }
      closeDialog(top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
 }

   initializeDialog(title=gettextRcmdr("Remove columns ..."))
    selFrame <- ttklabelframe(top, text=gettextRcmdr("Permanently remove column(s) from design"))
    estbuttonFrame <- ttkframe(selFrame)
    selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"),
            foreground = "darkgreen", command = onSelect,
            default = "normal", borderwidth = 3)
    tkgrid(selectButton)
    deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"),
            foreground = "darkgreen", command = onDeselect,
            default = "normal", borderwidth = 3)
    tkgrid(deselectButton)

    putRcmdr("sel.columns", variableListBox(selFrame, variableList=currem, listHeight=10, 
        title="Columns that will be removed",selectmode="multiple"))
    putRcmdr("other.removables", variableListBox(selFrame, variableList=removable, listHeight=10, 
        title="Columns that are removable",selectmode="multiple"))

         tkconfigure(sel.columns$listbox, listvariable=tclVar(paste(currem,collapse=" ")))
         sel.columns$varlist <- currem
         tkconfigure(other.removables$listbox, listvariable=tclVar(paste(removable,collapse=" ")))
         other.removables$varlist <- removable
         tkgrid(getFrame(other.removables), estbuttonFrame, getFrame(sel.columns))
         tkgrid(selFrame,sticky="n")

    OKCancelHelp(helpSubject="Menu.colremove")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=3)
}