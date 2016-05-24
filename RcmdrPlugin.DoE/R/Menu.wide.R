## one instance of assign replaced
Menu.wide <- function(){
initializeDialog(window=tab6, title=gettextRcmdr("Change from long to wide format ..."),
       offset=c(-commanderPosition()+10))   
       ## my offset is ignored if commander itself has upper left corner outside screen

onOK <- function(){
     closeDialog(window=tab6)
     .activeDataSet <- activeDataSet()
     di <- design.info(eval(parse(text=ActiveDataSet())))
  ### capture error messages from export function
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(window=tab6,recall=Menu.wide, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (length(grep("param",di$type))>0) 
            command <- paste(name," <- paramtowide(", .activeDataSet, ")",sep="")
        else command <- paste(name," <- reptowide(", .activeDataSet, ")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6,recall=Menu.wide, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(name, hilf, envir=.GlobalEnv)
        justDoIt(paste(name, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        activeDataSet(name)
        closeDialog(window=tab6)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

     namechange <- function(){
        if (!is.valid.name(tclvalue(nameVar)))
        tkmessageBox(message="invalid name!",icon="error", type="ok", 
        title="Invalid name for wide design")
    }
 ######## end define functions                          

##### define userform


## tab6 for exporting
#helptab6Button <- buttonRcmdr(tab6, text = gettextRcmdr("Tab Help"), 
#        foreground = "darkgreen", command = onHelpTab6, 
#        default = "normal", borderwidth = 3)

.activeDataSet <- ActiveDataSet()
putRcmdr("nameVar", tclVar(paste(.activeDataSet,"wide",sep=".")))

## settings for aggregation function
newFrame <- ttklabelframe(tab6, text=gettextRcmdr("Name of wide data frame:"))
newnamEntry <- tkentry(newFrame, width="50", textvariable=nameVar)
tkbind(newnamEntry, "<FocusOut>", namechange)
tkgrid(newnamEntry, sticky="w")
tkgrid(newFrame, sticky="w")

OKCancelHelp(window=tab6, helpSubject="Menu.exportTab")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=tab6, rows=5, columns=3, focus=newnamEntry, bindReturn=FALSE)

}
###
# End of Menu.export
###
