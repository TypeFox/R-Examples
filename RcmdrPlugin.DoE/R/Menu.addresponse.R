## one instance of assign replaced by justDoIt

Menu.addresponse <- function(){
initializeDialogDoE(title=gettextRcmdr("Add response ..."))   
       ## refresh bei Nutzung der radiobuttons
       ## Link auf das Menü des R-Commanders zum Einbinden von Daten --> anschl. R-Objekt verwenden

if (!exists("nameVar")) putRcmdr("nameVar", tclVar(""))
if (!exists("newnameVar")) putRcmdr("newnameVar", tclVar(""))
if (!exists("respVar")) putRcmdr("respVar", tclVar(""))
if (!exists("responseQuelle")) putRcmdr("responseQuelle", "R")

onOK <- function(){
     closeDialog(window=topdes2)
  ### capture error messages from export function
        name <- tclvalue(nameVar)
        newname <- tclvalue(newnameVar)
        putRcmdr("csvpath", tclvalue(fileVar))
        putRcmdr("respname", tclvalue(respVar))
        putRcmdr("decimal.setting", tclvalue(decimalrbVariable))
        putRcmdr("responseQuelle", tclvalue(etyperbVariable))
        if (!identical(respname,"")){
              hilf <- try(eval(parse(text=tclvalue(respVar))))
              if (class(hilf)[1]=="try-error"){ 
                  errorCondition(window=topdes2,recall=Menu.addresponse, 
                  gettextRcmdr("Invalid response specification"))
                  return()
                  }
              putRcmdr(respname, hilf)
              }
        if (is.element(newname, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(newname, gettextRcmdr("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.addresponse, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
        if (!is.valid.name(name)) {
            errorCondition(window=topdes2,recall=Menu.addresponse, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (!is.valid.name(newname)) {
            errorCondition(window=topdes2,recall=Menu.addresponse, 
                    message=paste('"', newname, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
    ### adding a response
        if (tclvalue(etyperbVariable)=="R")
           command <- paste("add.response(",name,
               ",", respname, ", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
           else{
              if (tclvalue(decimalrbVariable)=="default") command <- paste("add.response(",name,
                     ", ",dquote(csvpath),", replace=",
                     as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
              else command <- paste("add.response(",name, 
                     ", ",dquote(csvpath),", replace=",
                     as.logical(as.numeric(tclvalue(replacecbVariable))),", InDec=", 
                     dquote(tclvalue(decimalrbVariable)),")",sep="")
               }
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error"){
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
             errorCondition(window=topdes2, recall=Menu.addresponse, message=gettextRcmdr(hilf))
             return()
            }
        ## replace assign by justDoIt; assign(newname, hilf, envir=.GlobalEnv)
        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(newname, hilf, envir=.GlobalEnv)
        justDoIt(paste(newname, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        logger(paste(newname, "<-", command))
        activeDataSet(newname)
        putRcmdr("nameVar", tclVar(""))
        putRcmdr("newnameVar", tclVar(""))
        putRcmdr("csvpath", "")
        putRcmdr("respVar", tclVar(""))
        closeDialog(window=topdes2)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

     namechange <- function(){
        if (!tclvalue(respVar)==""){
        if (!exists(tclvalue(respVar)))
          tkmessageBox(message="invalid response name!", icon="error", type="ok", title="Non-existing response")}
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


 onRadio <- function(){
    if (tclvalue(etyperbVariable)=="R"){
        tkconfigure(dirEntry, state="disabled")
        tkconfigure(dirButton, state="disabled")
        tkconfigure(fileEntry, state="disabled")
        tkconfigure(fileButton, state="disabled")
        tkconfigure(defaultrb, state="disabled")
        tkconfigure(pointrb, state="disabled")
        tkconfigure(commarb, state="disabled")
        tkconfigure(respEntry, state="normal")
    } 
    else{ 
        tkconfigure(dirEntry, state="normal")
        tkconfigure(dirButton, state="normal")
        tkconfigure(fileEntry, state="normal")
        tkconfigure(fileButton, state="normal")
        tkconfigure(defaultrb, state="normal")
        tkconfigure(pointrb, state="normal")
        tkconfigure(commarb, state="normal")
        tkconfigure(respEntry, state="disabled")
    } 
 }
 
 onChangeDir <- function(){
     putRcmdr("direct",tclvalue(tkchooseDirectory()))
     if (!direct=="") {
        putRcmdr("dirVar", tclVar(direct))
        tkconfigure(dirEntry, textvariable = dirVar)
        setwd(direct)
        }
 }
 
 onChangeFile <- function(){
     fn <- tclvalue(tkgetOpenFile(filetypes=
        gettextRcmdr('{"csv files" {".csv" ".Csv" ".CSV"}} {"All Files" {"*"}}')))
     if (!fn=="") {
        putRcmdr("fileVar", tclVar(fn))
        tkconfigure(fileEntry, textvariable = fileVar)
        }
 }
 ######## end define functions                          

##### define userform


## topdes2 for adding a response
#helptopdes2Button <- buttonRcmdr(topdes2, text = gettextRcmdr("Tab Help"), 
#        foreground = "darkgreen", command = onHelptopdes2, 
#        default = "normal", borderwidth = 3)

putRcmdr("nameVar", tclVar(""))
putRcmdr("newnameVar", tclVar(""))
designs <- listDesigns()

designFrame <- tkframe(topdes2)
designsel <- ttkcombobox(designFrame, textvariable=nameVar, values=designs, state="readonly")
tcl(designsel, "current", 0)
.activeDataSet <- ActiveDataSet()

if (!is.null(.activeDataSet))
if (.activeDataSet %in% designs) tcl(designsel, "current", which(designs==.activeDataSet)-1)
    tkbind(designsel, "<FocusOut>", namechange)

replacecbVariable <- tclVar("0")
replacecb <- ttkcheckbutton(designFrame,text=gettextRcmdr("Replace responses, if they exist already"),
      variable=replacecbVariable)

tkgrid(ttklabel(designFrame,text="Choose design to which one or more responses are to be added:"),sticky="w",columnspan=3) 
tkgrid(designsel, replacecb, sticky="w", padx=15)

#putRcmdr("nameVar", tclVar(designs[as.numeric(tclvalue(tcl(designsel, "current")))+1]))

## radio buttons for choosing response source type
etradioFrame <- ttklabelframe(topdes2, text=gettextRcmdr("What type of response ?"))
etyperbVariable <- tclVar(responseQuelle)
Rrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("R object (data frame or vector"),variable=etyperbVariable,value="R",command=onRadio)
csvrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("External csv file"),variable=etyperbVariable,value="csv",command=onRadio)
tkgrid(Rrb, sticky="w")
tkgrid(csvrb, sticky="w")


## import file
## (enabled visible for external only)
fileFrame <- ttklabelframe(topdes2, text=gettextRcmdr("csv file with response data"))
putRcmdr("dirVar", tclVar(getwd()))
dirEntry <- tkentry(fileFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(fileFrame, text = gettextRcmdr("Change \n working directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(fileFrame, text="   "), dirButton, sticky="w")
## radio buttons for choosing export decimal separator
## make visible for external only!!!
decimalradioFrame <- ttklabelframe(fileFrame, text=gettextRcmdr("Decimal Separator ?"))
decimalrbVariable <- tclVar("default")
if (exists("decimal.setting")) decimalrbVariable <- tclVar(decimal.setting)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

if (!exists("csvpath")) putRcmdr("fileVar", tclVar(""))
   else putRcmdr("fileVar", tclVar(csvpath))
fileEntry <- tkentry(fileFrame, width="50", textvariable=fileVar)
fileButton <- buttonRcmdr(fileFrame, text = gettextRcmdr("Select csv file"), 
        foreground = "darkgreen", width = "20", command = onChangeFile, 
        default = "normal", borderwidth = 3)
tkgrid(fileEntry, tklabel(fileFrame, text="   "), fileButton, sticky="w")
tkgrid(decimalradioFrame, sticky="w", padx=50, pady=10)

## response object within R
## can be numeric vector or matrix or data frame
##    numeric vector can be calculated on the fly
respFrame <- tkframe(topdes2)
respEntry <- tkentry(respFrame, width="20", textvariable=respVar)
efnamelabel <- tklabel(respFrame,text=gettextRcmdr("R object that contains the response(s)"))
tkgrid(efnamelabel, respEntry, sticky="w", pady=20)
tkgrid.configure(respEntry, padx=15)

newnamelabel <- tklabel(topdes2, text="Name for new design")
putRcmdr("newnameVar", tclVar(tclvalue(nameVar)))
newnameEntry <- tkentry(topdes2, width="50", textvariable=newnameVar)

tkgrid(designFrame, sticky="w",pady=15, columnspan=3)

tkgrid(etradioFrame, sticky="w",pady=15)

tkgrid(respFrame, columnspan=3, sticky="w")

tkgrid(fileFrame, sticky="w", columnspan=3, pady=15)


## default: R object
#        tkconfigure(fileEntry, state="disabled")
#        tkconfigure(fileButton, state="disabled")
#        tkconfigure(defaultrb, state="disabled")
#        tkconfigure(pointrb, state="disabled")
#        tkconfigure(commarb, state="disabled")
onRadio()

tkgrid(newnamelabel, sticky="w")
tkgrid(newnameEntry, sticky="w")

OKCancelHelp(window=topdes2, helpSubject="Menu.addresponse")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=topdes2, rows=5, columns=3, focus=designsel, bindReturn=FALSE)

}
###
# End of Menu.addresponse
###
