## one instance of assign replaced by justDoIt

Menu.importrdacsv <- function(){
## Aktualisierung des neuen Namens bei Änderung in designsel funktioniert noch nicht
## tkbind event vermutlich nicht gut gewählt!

if (!exists("nameVar")) putRcmdr("nameVar", tclVar(""))
if (!exists("newnameVar")) putRcmdr("newnameVar", tclVar(""))
if (!exists("loadVar")) putRcmdr("designs", character(0))
   else {
   if (tclvalue(loadVar)=="") putRcmdr("designs", character(0)) 
   else putRcmdr("designs", listDesigns())
   }

onOK <- function(){
     closeDialog(window=topdes2)
  ### capture error messages from export function
        name <- tclvalue(nameVar)
        newname <- tclvalue(newnameVar)
        putRcmdr("decimal.setting", tclvalue(decimalrbVariable))
        if (is.element(newname, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(newname, gettextRcmdr("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.importrdacsv, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
        if (!is.valid.name(name)) {
            errorCondition(window=topdes2,recall=Menu.importrdacsv, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (!is.valid.name(newname)) {
            errorCondition(window=topdes2,recall=Menu.importrdacsv, 
                    message=paste('"', newname, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
    ### adding a response
        putRcmdr("csvpath", tclvalue(fileVar))
        name <- tclvalue(nameVar)
        newname <- tclvalue(newnameVar)
              if (tclvalue(decimalrbVariable)=="default") command <- paste("add.response(",name,
                     ", ",dquote(csvpath), ", replace=",
                     as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
              else command <- paste("add.response(",name,
                     ", ",dquote(csvpath), ", replace=",
                     as.logical(as.numeric(tclvalue(replacecbVariable))),", InDec=",
                     dquote(tclvalue(decimalrbVariable)),")",sep="")

      hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.importrdacsv, message=gettextRcmdr(hilf))
             return()
            }

        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(newname, hilf, envir=.GlobalEnv)
        justDoIt(paste(newname, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        logger(paste(newname, "<-", command))
        activeDataSet(newname)
        putRcmdr("fileVar", tclVar(""))
        putRcmdr("loadVar", tclVar(""))
        putRcmdr("nameVar", tclVar(""))
        putRcmdr("newnameVar", tclVar(""))
        closeDialog(window=topdes2)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
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

 onChangerdaFile <- function(){
     fn <- tclvalue(tkgetOpenFile(filetypes=
        gettextRcmdr('{"R Data Files" {".rda" ".Rda" ".RDA" ".RData"}} {"All Files" {"*"}}')))
     if (!fn=="") {
        putRcmdr("loadVar", tclVar(fn))
        tkconfigure(loadEntry, textvariable = loadVar)
        command <- paste("load(",dQuote(tclvalue(loadVar)),")")
        justDoItDoE(command)
        logger(command)
        putRcmdr("designs", listDesigns())
        tkconfigure(designsel, values=designs)
        if (length(designs)>0) tcl(designsel, "current", 0)
        }
 }

 onChangecsvFile <- function(){
     fn <- tclvalue(tkgetOpenFile(filetypes=
        gettextRcmdr('{"csv files" {".csv" ".Csv" ".CSV"}} {"All Files" {"*"}}')))
     if (!fn=="") {
        putRcmdr("fileVar", tclVar(fn))
        tkconfigure(fileEntry, textvariable = fileVar)
        }
 }
     namechange <- function(){
        if (tclvalue(newnameVar)=="" & !tclvalue(nameVar)=="")
        putRcmdr("newnameVar", paste(tclvalue(nameVar),"withresp",sep="."))
        tkconfigure(newnameEntry, textvariable=newnameVar)
    }
 ######## end define functions                          

##### define userform
initializeDialogDoE(title=gettextRcmdr("Import design from rda and csv ..."))   

dirFrame <- ttklabelframe(topdes2, text=gettextRcmdr("Directory"))
putRcmdr("dirVar", tclVar(getwd()))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change \n working directory"), 
        foreground = "darkgreen", width = "20", command=onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")
tkgrid(dirFrame)

## import file names for rda and csv
loadFrame <- ttklabelframe(topdes2, text=gettextRcmdr("rda file with design"))
if (!exists("loadVar")) putRcmdr("loadVar", tclVar(""))
loadEntry <- tkentry(loadFrame, width="50", textvariable=loadVar)
loadButton <- buttonRcmdr(loadFrame, text = gettextRcmdr("Select rda file"), 
        foreground = "darkgreen", width = "20", command = onChangerdaFile, 
        default = "normal", borderwidth = 3)
tkgrid(loadEntry, tklabel(loadFrame, text="   "), loadButton, sticky="w")


fileFrame <- ttklabelframe(topdes2, text=gettextRcmdr("csv file with response data"))
if (!exists("fileVar")) putRcmdr("fileVar", tclVar(""))
fileEntry <- tkentry(fileFrame, width="50", textvariable=fileVar)
fileButton <- buttonRcmdr(fileFrame, text = gettextRcmdr("Select csv file"), 
        foreground = "darkgreen", width = "20", command = onChangecsvFile, 
        default = "normal", borderwidth = 3)
      ## radio buttons for choosing import decimal separator
decimalradioFrame <- ttklabelframe(fileFrame, text=gettextRcmdr("Decimal Separator ?"))
decimalrbVariable <- tclVar("default")
if (exists("decimal.setting")) decimalrbVariable <- tclVar(decimal.setting)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")
tkgrid(fileEntry, fileButton, sticky="w")
tkgrid.configure(fileButton, padx=5)
tkgrid(decimalradioFrame, sticky="w", padx=50, pady=10)


## old design name
designFrame <- tkframe(topdes2)
designsel <- ttkcombobox(designFrame, textvariable=nameVar, values=designs, state="readonly")
if (length(designs)==1) tcl(designsel, "current", 0)
.activeDataSet <- ActiveDataSet()
    tkbind(designsel, "<<ComboboxSelected>>", namechange)

replacecbVariable <- tclVar("0")
replacecb <- ttkcheckbutton(designFrame,text=gettextRcmdr("Replace responses, \nif they exist already"),
      variable=replacecbVariable)

if (!is.null(.activeDataSet))
if (.activeDataSet %in% designs) tcl(designsel, "current", which(designs==.activeDataSet)-1)
tkgrid(designsel, replacecb, sticky="w")

newnamelabel <- tklabel(topdes2, text="Name for new design")
if (tclvalue(newnameVar)=="" & !tclvalue(nameVar)=="") 
     putRcmdr("newnameVar", tclVar(paste(tclvalue(nameVar),"withresp",sep=".")))
newnameEntry <- tkentry(topdes2, width="50", textvariable=newnameVar)


## always grid details, as otherwise default file name does not work
## design name info and help button have already been gridded above
tkgrid(loadFrame, sticky="w", columnspan=3, pady=20)
tkgrid(tklabel(topdes2, text=gettextRcmdr("Select design")), sticky="w", columnspan=3)
tkgrid(designFrame, sticky="w", columnspan=3)
tkgrid.configure(replacecb, padx=15)
tkgrid(fileFrame, sticky="w", columnspan=3, pady=20)
tkgrid(newnamelabel, sticky="w")
tkgrid(newnameEntry, sticky="w")

OKCancelHelp(window=topdes2, helpSubject="Menu.importrdacsv")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="s", columnspan=3)

dialogSuffix(window=topdes2, rows=5, columns=3, focus=newnameEntry, bindReturn=FALSE)

}
###
# End of Menu.addresponse
###
