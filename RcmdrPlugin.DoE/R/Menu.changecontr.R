Menu.changecontr <- function(){
    ## view and change response names
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))

    putRcmdr("curch", NULL)
    quant <- di$quantitative
    if (!is.null(quant)) 
    putRcmdr("changeable", names(di$factor.names)[!quant])
    else 
    putRcmdr("changeable", names(di$factor.names))

   onSelect <- function(){
     if (tclvalue(tcl(other.changeables$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- changeable[as.numeric(unlist(strsplit(tclvalue(tcl(other.changeables$listbox, "curselection")), " ")))+1]
     putRcmdr("changeable", setdiff(changeable,add))
     putRcmdr("curch", c(getRcmdr("curch"),add))
     add <- NULL
     tcl(other.changeables$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.changeables$listbox, listvariable=tclVar(paste(changeable,collapse=" ")))
     other.changeables$varlist <- changeable
     tkconfigure(sel.tochange$listbox, listvariable=tclVar(paste(curch,collapse=" ")))
     sel.tochange$varlist <- curch
   }
   onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(sel.tochange$listbox, "curselection"))=="") return()
     add <- curch[as.numeric(unlist(strsplit(tclvalue(tcl(sel.tochange$listbox, "curselection")), " ")))+1]
     putRcmdr("curch", setdiff(curch,add))
     putRcmdr("changeable", c(getRcmdr("changeable"),add))
     add <- NULL
     tcl(sel.tochange$listbox, "selection", "clear", "0", "999")
     tkconfigure(other.changeables$listbox, listvariable=tclVar(paste(changeable,collapse=" ")))
     other.changeables$varlist <- changeable
     tkconfigure(sel.tochange$listbox, listvariable=tclVar(paste(curch,collapse=" ")))
     sel.tochange$varlist <- curch
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
      hilf <- curch
      if (length(hilf)>0){ 
        command <- paste(ActiveDataSet(), "<- change.contr(",ActiveDataSet(),", c(", 
             paste(hilf,rep(dquote(tclvalue(contr_rbVariable)),length(hilf)),sep="=",collapse=","),"))")
        doItAndPrint(command)
        putRcmdr(".activeDataSet", ActiveDataSet())
        activeDataSet(.activeDataSet)
      }
      closeDialog(top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
 }

   initializeDialog(title=gettextRcmdr("Change contrasts ..."))
    selFrame <- ttklabelframe(top, text=gettextRcmdr("Change contrasts for design factors ..."))
    estbuttonFrame <- ttkframe(selFrame)
    selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"),
            foreground = "darkgreen", command = onSelect,
            default = "normal", borderwidth = 3)
    tkgrid(selectButton)
    deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"),
            foreground = "darkgreen", command = onDeselect,
            default = "normal", borderwidth = 3)
    tkgrid(deselectButton)

    putRcmdr("sel.tochange", variableListBox(selFrame, variableList=curch, listHeight=10, 
        title="Contrasts will be changed",selectmode="multiple"))
    putRcmdr("other.changeables", variableListBox(selFrame, variableList=changeable, listHeight=10, 
        title="Contrasts will not be changed",selectmode="multiple"))
        
         tkconfigure(sel.tochange$listbox, listvariable=tclVar(paste(curch,collapse=" ")))
         sel.tochange$varlist <- curch
         tkconfigure(other.changeables$listbox, listvariable=tclVar(paste(changeable,collapse=" ")))
         other.changeables$varlist <- changeable
         tkgrid(getFrame(other.changeables), estbuttonFrame, getFrame(sel.tochange))
         tkgrid(selFrame,sticky="n")
         

    contrFrame <- ttklabelframe(selFrame, text=gettextRcmdr("Contrast type ..."))
putRcmdr("contr_rbVariable", tclVar("contr.treatment"))
trt_rb <- tkradiobutton(contrFrame,text=gettextRcmdr("Treatment (dummy) contrasts"),variable=contr_rbVariable,value="contr.treatment")
FrF2_rb <- tkradiobutton(contrFrame,text=gettextRcmdr("FrF2 contrasts (number of levels must be power of 2)"),variable=contr_rbVariable,value="contr.FrF2")
sum_rb <- tkradiobutton(contrFrame,text=gettextRcmdr("Helmert contrasts"),variable=contr_rbVariable,value="contr.helmert")
helmert_rb <- tkradiobutton(contrFrame,text=gettextRcmdr("Sum (deviation) contrasts"),variable=contr_rbVariable,value="contr.sum")
poly_rb <- tkradiobutton(contrFrame,text=gettextRcmdr("Polynomial contrasts"),variable=contr_rbVariable,value="contr.poly")

tkgrid(trt_rb, sticky="w")
tkgrid(FrF2_rb, sticky="w")
tkgrid(sum_rb, sticky="w")
tkgrid(helmert_rb, sticky="w")
tkgrid(poly_rb, sticky="w")

tkgrid(contrFrame, sticky="w", columnspan=3)

    OKCancelHelp(helpSubject="Menu.changecontr")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=3)
}