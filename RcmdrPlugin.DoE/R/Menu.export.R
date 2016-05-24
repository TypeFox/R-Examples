Menu.export <- function(){
initializeDialogDoE(title=gettextRcmdr("Export design ..."))   
       ## my offset is ignored if commander itself has upper left corner outside screen

onOK <- function(){
     closeDialog(window=topdes2)
  ### capture error messages from export function
        name <- tclvalue(nameVar)
     putRcmdr("decimal.setting", tclvalue(decimalrbVariable))
        if (!is.valid.name(name)) {
            errorCondition(window=topdes2,recall=Menu.export, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
    ### exporting
        putRcmdr("path", tclvalue(dirVar))
        putRcmdr("filename", tclvalue(fileVar))
        if (tclvalue(decimalrbVariable)=="default") command <- paste("export.design(",name,
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(getRcmdr("path")),", file=",dquote(getRcmdr("filename")),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
        else command <- paste("export.design(",name, 
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(getRcmdr("path")),", file=",dquote(getRcmdr("filename")),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),", OutDec=", 
               dquote(tclvalue(decimalrbVariable)),")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.export, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        activeDataSet(name)
        closeDialog(window=topdes2)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

     nameenter <- function(){
           if (identical(tclvalue(getRcmdr("fileVar")),tclvalue(getRcmdr("nameVar")))
              || tclvalue(getRcmdr("fileVar"))=="")
              putRcmdr("name.equal.filename", TRUE)
           else putRcmdr("name.equal.filename", FALSE)
        }
     namechange <- function(){
        if (is.valid.name(tclvalue(nameVar))){
          if (name.equal.filename){
          putRcmdr("fileVar", tclVar(tclvalue(nameVar)))  ## otherwise, variables would be directly tied
          tkconfigure(fileEntry, textvariable=getRcmdr("fileVar"))
          }
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid design name")
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

 onChangeDir <- function(){
     putRcmdr("direct",tclvalue(tkchooseDirectory()))
     if (!direct=="") {
        putRcmdr("dirVar", tclVar(direct))
        tkconfigure(dirEntry, textvariable = dirVar)
        setwd(direct)
        }
 }
 ######## end define functions                          

##### define userform


## topdes2 for exporting
#helptopdes2Button <- buttonRcmdr(topdes2, text = gettextRcmdr("Tab Help"), 
#        foreground = "darkgreen", command = onHelptopdes2, 
#        default = "normal", borderwidth = 3)

putRcmdr("nameVar", tclVar(""))
designs <- listDesigns()
designsel <- ttkcombobox(topdes2, textvariable=nameVar, values=designs, state="readonly")
tcl(designsel, "current", 0)
.activeDataSet <- ActiveDataSet()

if (!is.null(.activeDataSet))
if (.activeDataSet %in% designs) tcl(designsel, "current", which(designs==.activeDataSet)-1)
    tkbind(designsel, "<FocusIn>", nameenter)
    tkbind(designsel, "<FocusOut>", namechange)

#putRcmdr("nameVar", tclVar(designs[as.numeric(tclvalue(tcl(designsel, "current")))+1]))

tkgrid(ttklabel(topdes2,text="Choose design to be saved:"),designsel,sticky="w") 
tkgrid(ttklabel(topdes2,text="   "),sticky="w") 

## radio buttons for choosing export type
etradioFrame <- ttklabelframe(topdes2, text=gettextRcmdr("(How to) Export ?"))
etyperbVariable <- tclVar("html")
allrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("all formats"),variable=etyperbVariable,value="all")
rdarb <- tkradiobutton(etradioFrame,text=gettextRcmdr("rda only"),variable=etyperbVariable,value="rda")
htmlrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("html and rda"),variable=etyperbVariable,value="html")
csvrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("csv and rda"),variable=etyperbVariable,value="csv")
tkgrid(allrb, sticky="w")
tkgrid(rdarb, sticky="w")
tkgrid(htmlrb, sticky="w")
tkgrid(csvrb, sticky="w")

## radio buttons for choosing export decimal separator
decimalradioFrame <- ttklabelframe(topdes2, text=gettextRcmdr("Decimal Separator ?"))
decimalrbVariable <- tclVar("default")
if (exists("decimal.setting", mode="character")) decimalrbVariable <- tclVar(decimal.setting)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(topdes2, text=gettextRcmdr("Storage Directory"))
putRcmdr("dirVar", tclVar(getwd()))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change \n working directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")

## export file name
putRcmdr("fileVar", tclVar(tclvalue(nameVar)))
fileEntry <- tkentry(topdes2, width="20", textvariable=fileVar)
efnamelabel <- tklabel(topdes2,text=gettextRcmdr("Export file names: name below with appropriate endings (html or csv, and rda)"))
replacecbVariable <- tclVar("0")
replacecb <- ttkcheckbutton(topdes2,text=gettextRcmdr("Replace file(s), if exists"),variable=replacecbVariable)

## always grid details, as otherwise default file name does not work
## design name info and help button have already been gridded above
tkgrid(etradioFrame, decimalradioFrame, sticky="w")
tkgrid(tklabel(topdes2,text="   "))
tkgrid(dirFrame, sticky="w", columnspan=5)
tkgrid(tklabel(topdes2,text="  "), sticky="w")
tkgrid(efnamelabel, sticky="w", columnspan=5)
tkgrid(fileEntry, sticky="w", columnspan=5)
tkgrid(replacecb, sticky="w", columnspan=5)

OKCancelHelp(window=topdes2, helpSubject="Menu.exportTab")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=topdes2, rows=5, columns=3, focus=designsel, bindReturn=FALSE)

}
###
# End of Menu.export
###
