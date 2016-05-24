Menu.pickcube <- function(){
initializeDialog(window=tab6, title=gettextRcmdr("Pick cube portion from central composite design ..."),
       offset=c(-commanderPosition()+10))   
       ## my offset is ignored if commander itself has upper left corner outside screen

onOK <- function(){
     closeDialog(window=tab6)
  ### capture error messages from export function
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(window=tab6,recall=Menu.pickcube, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        command <- paste(tclvalue(nameVar), " <- pickcube(",.activeDataSet, ")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6,recall=Menu.pickcube, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        activeDataSet(name)
        closeDialog(window=tab6)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

     namechange <- function(){
        if (!is.valid.name(tclvalue(nameVar)))
        tkmessageBox(message="invalid name!",icon="error", type="ok", 
        title="Invalid name for reduced (cube only) design")
    }

 ######## end define functions                          

##### define userform


## tab6 for exporting
#helptab6Button <- buttonRcmdr(tab6, text = gettextRcmdr("Tab Help"), 
#        foreground = "darkgreen", command = onHelpTab6, 
#        default = "normal", borderwidth = 3)

.activeDataSet <- ActiveDataSet()

## settings for aggregation function
putRcmdr("nameVar", tclVar(paste(.activeDataSet, "cubeonly", sep=".")))
newFrame <- ttklabelframe(tab6, text=gettextRcmdr("Name of reduced data frame:"))
newnamEntry <- tkentry(newFrame, width="50", textvariable=nameVar)

tkgrid(tklabel(tab6, text=gettextRcmdr("Name of active design:")), tklabel(tab6, text=.activeDataSet), sticky="w")

tkbind(newnamEntry, "<FocusOut>", namechange)
tkgrid(newnamEntry, sticky="w")
tkgrid(newFrame, sticky="w")

OKCancelHelp(window=tab6, helpSubject="Menu.pickcube")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=tab6, rows=5, columns=3, focus=newnamEntry, bindReturn=FALSE)

}
###
# End of Menu.pickcube
###
