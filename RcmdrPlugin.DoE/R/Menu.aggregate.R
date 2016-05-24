Menu.aggregate <- function(){
initializeDialog(window=tab6, title=gettextRcmdr("Change from long to wide format ..."),
       offset=c(-commanderPosition()+10))   
       ## my offset is ignored if commander itself has upper left corner outside screen

onOK <- function(){
     closeDialog(window=tab6)
  ### capture error messages from export function
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(window=tab6,recall=Menu.aggregate, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
   ##     if (tclvalue(transfoVar)=="none") transfoVar <- tclVar("NULL")
   ##     command <- paste(tclvalue(nameVar), " <- aggregate(",.activeDataSet, 
   ##            ", response=",dquote(tclvalue(respVar)),",transformation=",dquote(tclvalue(transfoVar)),", FUN=",dquote(tclvalue(funVar)),")",sep="")
        command <- paste(tclvalue(nameVar), " <- aggregate(",.activeDataSet, 
               ", response=",dquote(tclvalue(respVar)),", FUN=",dquote(tclvalue(funVar)),")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6,recall=Menu.aggregate, message=gettextRcmdr(hilf))
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
        title="Invalid name for aggregate design")
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

 ######## end define functions                          

##### define userform


## tab6 for exporting
#helptab6Button <- buttonRcmdr(tab6, text = gettextRcmdr("Tab Help"), 
#        foreground = "darkgreen", command = onHelpTab6, 
#        default = "normal", borderwidth = 3)

.activeDataSet <- ActiveDataSet()

## settings for aggregation function
putRcmdr("nameVar", tclVar(.activeDataSet))
newFrame <- ttklabelframe(tab6, text=gettextRcmdr("Name of aggregated data frame:"))
newnamEntry <- tkentry(newFrame, width="50", textvariable=nameVar)
tkbind(newnamEntry, "<FocusOut>", namechange)
tkgrid(newnamEntry, sticky="w")
tkgrid(newFrame, sticky="w")
funVar <- tclVar("mean")
FUNEntry <- ttkcombobox(tab6, textvariable=funVar, values=c("mean","sd","SN"), state="normal")
tkgrid(tklabel(tab6, text="Function used for aggregation:"), sticky="w")
tkgrid(FUNEntry, sticky="w")
respVar <- tclVar(colnames(design.info(eval(parse(text=.activeDataSet)))$responselist)[1])
respEntry <- ttkcombobox(tab6, textvariable=respVar, 
     values=colnames(design.info(eval(parse(text=.activeDataSet)))$responselist), state="readonly")
tkgrid(tklabel(tab6, text="Response to be aggregated:"), sticky="w")
tkgrid(respEntry, sticky="w")
##transfoVar <- tclVar("none")
##transfoEntry <- tkentry(tab6, width="50", textvariable=transfoVar)
##tkgrid(tklabel(tab6, text="Transformation of raw data before aggregation, if requested:"), sticky="w")
##tkgrid(transfoEntry, sticky="w")


OKCancelHelp(window=tab6, helpSubject="Menu.aggregate")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=tab6, rows=5, columns=3, focus=newnamEntry, bindReturn=FALSE)

}
###
# End of Menu.export
###
