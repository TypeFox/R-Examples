## one instance of assign replaced by justDoItDoE with modified command and justDoIt("design.info(name) <- hilfatt")
## one instance of assign replaced by justDoIt

Menu.Dopt <- function(){
    putRcmdr(".activeDataSet", ActiveDataSet())
    putRcmdr("variables", colnames(get(getRcmdr(".activeDataSet"))))
#candidate set is the active data set

initializeDialogDoE(title=gettextRcmdr("D-optimal design ..."))   
     ## function initializeDialogDoE assumes topdes2 as windowname
     ## last stored top left corner for window is stored under topleft2xy
     ## onRefresh still makes window walk a little

if (exists("curindex", where="RcmdrEnv")) rm(curindex, pos="RcmdrEnv")

if (!exists(".stored.designDopt",where="RcmdrEnv")) 
           putRcmdr(".stored.designDopt", .default.designDopt)
           ## nameVar, nrunVar, nfacVar, nrepVar
           ## cbInitials containing repeat.onlyVariable, randomizeVariable, 
           ##                       facnamesAutoVariable, faclevelsCommonVariable, 
           ##                       nrunEntryVariable, estcbVariable
           ##                       specialcbVariable, replacecbVariable, MaxC2cbVariable
           ##                       res3cbVariable
           ## level1Var, level2Var, seedVar, specialrbVariable, hardVar, genVar, 
           ## catlgVar, designVar, designrbVariable, destyperbVariable
           ## resVar, qualcritrbVariable, facnamlist,faclev1list,faclev2list, faclablist
           ## etyperbVariable, decimalrbVariable, dirVar, fileVar

## MaxC2cbVariable is free again (no. 9 of cbInitials)

## define called functions
 infoClose <- function(){
     putRcmdr("infotxt",tclVar(""))
 }
 
 onHelpTab1 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.DoptTab1"))
 }
 onHelpTab6 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.exportTab"))
 }
 
 tabpos <- function(){
        ### get 0-based index of currently selected tab
        activestab.tn <- tclvalue(tcl(tn, "select"))
        activestab.tn <- strsplit(activestab.tn,".",fixed=TRUE)[[1]]
        activestab.tn <- as.numeric(activestab.tn[length(activestab.tn)])-1
        activestab.tn
}

storeRcmdr <- function(){
        hilf <- list(nameVar=tclvalue(nameVar),
        nrunsVar=tclvalue(nrunsVar),
        rhsVariable=tclvalue(rhsVariable),
        constraintVar=tclvalue(constraintVar),
        nrepeatVar=tclvalue(nrepeatVar),
        cbInitials = c("0", tclvalue(randomizecbVariable),
                       "0",
                       1,0,
                       0,tclvalue(replacecbVariable),0,
                       0
                       ),
        seedVar=tclvalue(seedVar),
        etyperbVariable=tclvalue(etyperbVariable),
        decimalrbVariable=tclvalue(decimalrbVariable),
        dirVar=tclvalue(dirVar), fileVar=tclvalue(fileVar))
        class(hilf) <- c("menu.designDopt","list")
        putRcmdr(".stored.designDopt",hilf)
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
    onRefreshEnd()
    ## store entries so that users do not have to redo everything
    ## in case of stupid mistakes
    storeRcmdr()
    ## seed is not used from previously stored design
     closeDialog(window=topdes2)
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)){
            errorCondition(window=topdes2,recall=Menu.Dopt, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.Dopt, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
    ###  further error messages with return to menu ?

    ### not yet perfect, especially NULL entries are not possible
    ### for didactic reasons distinguish between usage of default.levels and other?
    if (!tclvalue(constraintVar) == gettextRcmdr("<all candidate data set rows eligible>"))
    command <- paste("Dopt.design(",tclvalue(nrunsVar), ", data=", getRcmdr(".activeDataSet"),
                  ", formula=", dquote(paste("~",tclvalue(rhsVariable),sep="")),", constraint=", dquote(tclvalue(constraintVar)),
                  ", nRepeat=",tclvalue(nrepeatVar), 
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizecbVariable))),
                  ",seed=",tclvalue(seedVar),")") 
    else     command <- paste("Dopt.design(",tclvalue(nrunsVar), ", data=", getRcmdr(".activeDataSet"),
                  ", formula=", dquote(paste("~",tclvalue(rhsVariable),sep="")),
                  ", nRepeat=",tclvalue(nrepeatVar), 
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizecbVariable))),
                  ",seed=",tclvalue(seedVar),")") 

        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.Dopt, message=gettextRcmdr(hilf))
             return()
            }
                  
        logger(paste(name, "<-", command))
        logger("## creator element of design.info will be different, when using the command line command!")
        ## change creator to contain menu settings
        hilfatt <- design.info(hilf)
        hilfatt$creator <- .stored.designDopt
        class(hilfatt$creator) <- c("menu.designDopt", "list")
        ## replace assign by justDoIt; assign(name, hilf, envir=.GlobalEnv)
        design.info(hilf) <- hilfatt
        putRcmdr("hilf", hilf)
        justDoIt(paste(name, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        activeDataSet(name)
    ### exporting
    if (!tclvalue(etyperbVariable)=="none"){
        putRcmdr("path", tclvalue(dirVar))
        putRcmdr("filename", tclvalue(fileVar))
        if (!as.logical(as.numeric(tclvalue(replacecbVariable)))){
          lf <- tolower(list.files(path = path))
          if (tolower(paste(filename, "rda", sep = ".")) %in% lf) 
                stop("file ", paste(filename, "rda", "."), 
                " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
          if (tclvalue(etyperbVariable)=="html" & tolower(paste(filename, "html", sep = ".")) %in% lf) 
                stop("file ", paste(filename, "html", "."), 
                " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
          if (tclvalue(etyperbVariable)=="csv" & tolower(paste(filename, "csv", sep = ".")) %in% lf) 
                stop("file ", paste(filename, "csv", "."), 
                " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
         }
        if (tclvalue(decimalrbVariable)=="default") command <- paste("export.design(",name,
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(path),", file=",dquote(filename),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
        else command <- paste("export.design(",name, 
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(path),", file=",dquote(filename),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),", OutDec=", dquote(tclvalue(decimalrbVariable)),")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            errorCondition(window=topdes2,recall=Menu.Dopt, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        }
        rm(activestab.tn, pos="RcmdrEnv")
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }


onLoad <- function(){
    ## seems to work now, needs to be tested!
        hilf <- listDesignlhs()
        if (length(hilf)==0) {
            tkmessageBox(message=gettextRcmdr("There are no stored design inputs in this session."),
            icon="error", type="ok", title="no stored design inputs")
            return()
            }
    putRcmdr("deschoose2",tktoplevel())
    tkwm.title(deschoose2, gettextRcmdr("Choose stored design form"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(deschoose2, position)
    putRcmdr("lb", variableListBox(deschoose2, variableList=hilf, title="Choose stored design form"))
        tkgrid(lb$frame)
    onOK <- function() {
        putRcmdr(".stored.designDopt",get(lb$varlist[as.numeric(tclvalue(tcl(lb$listbox, "curselection")))+1]))
        if ("design" %in% class(getRcmdr(".stored.designDopt"))) 
            putRcmdr(".stored.designDopt", design.info(getRcmdr(".stored.designDopt"))$creator)
        tkfocus(CommanderWindow())
        tkdestroy(topdes2)
        tkdestroy(deschoose2)
        Menu.Dopt()
    }
    OKCancelHelp(window=deschoose2)
    tkgrid(buttonsFrame, sticky="s")
    dialogSuffix(window=deschoose2, rows=1, columns=1, 
         focus=lb$listbox)
}

onRefreshEnd <- function(){
       storeRcmdr()
        ## letzte Position enthaelt tab index (beginnend bei 1)
        putRcmdr("activestab.tn",tabpos())
        ID <- topdes2$ID
        putRcmdr("topleft2xy",as.numeric(c(tclvalue(.Tcl(paste("winfo rootx", ID))), 
                              tclvalue(.Tcl(paste("winfo rooty", ID))))))
#        assign("activestab.tn",strsplit(activestab.tn,".",fixed=TRUE)[[1]],pos="RcmdrEnv")
#        assign("activestab.tn",as.numeric(activestab.tn[length(activestab.tn)])-1,pos="RcmdrEnv")
}

onRefresh <- function(){
#print(as.character(tclObj(tcl(tn, "select"))))
        onRefreshEnd()
        ## letzte Position enthaelt tab index (beginnend bei 1)
          tkfocus(CommanderWindow())
          tkdestroy(topdes2)
          Menu.Dopt()
}

onStore <- function(){
        ## Speichernamen abfragen und hier ermöglichen (statt stored.designDopt)
        textentry() ## creates text string stored in savename.RcmdrPlugin.DoE
        if (!is.null(savename.RcmdrPlugin.DoE)){
        if (!is.valid.name(savename.RcmdrPlugin.DoE)){
            textcorrect(gettextRcmdr("This is not a valid name. Please correct:"))
            return()
          }
        if (is.element(savename.RcmdrPlugin.DoE, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(savename.RcmdrPlugin.DoE, gettextRcmdr("Object"))))
            {
              textcorrect(gettextRcmdr("Please enter a new name:"))
              return()
             }
          }
        storeRcmdr()
        ## replace assign by justDoIt; assign(savename.RcmdrPlugin.DoE, getRcmdr(".stored.designDopt"), envir=.GlobalEnv)
        justDoIt(paste(savename.RcmdrPlugin.DoE, "<- getRcmdr(\".stored.designDopt\")"))
        message(gettextRcmdr("inputs have been stored"))
        }
}

onReset <- function(){
        assign(".stored.designDopt",.default.designDopt,pos="RcmdrEnv")
        tkfocus(CommanderWindow())
  tkdestroy(topdes2)
  Menu.Dopt()
}

     nameenter <- function(){
           if (identical(tclvalue(getRcmdr("fileVar")),tclvalue(getRcmdr("nameVar"))))
              putRcmdr("name.equal.filename", TRUE)
           else putRcmdr("name.equal.filename", FALSE)
        }
     namechange <- function(){
        if (is.valid.name(tclvalue(nameVar))){
          if (name.equal.filename){
          putRcmdr("fileVar", tclVar(tclvalue(nameVar)))  ## otherwise, variables would be directly tied
#          putRcmdr("exportlabVar", tclVar(paste("Current design to be saved:", tclvalue(nameVar),"\n   ")))
          tkconfigure(fileEntry, textvariable=getRcmdr("fileVar"))
#          tkconfigure(exportlab, textvariable=getRcmdr("exportlabVar"))
          }
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid design name")
    }


 onChangeDir <- function(){
     putRcmdr("direct",tclvalue(tkchooseDirectory()))
     if (!direct=="") {
        putRcmdr("dirVar", tclVar(direct))
        tkconfigure(dirEntry, textvariable = dirVar)
        }
 }


######## end define functions                          


##### define userform
#tn <- ttknotebook(top,height=100, width=500)


putRcmdr("tn",ttknotebook(topdes2))
#tn <- ttknotebook(topdes2)

putRcmdr("tab1",ttkframe(tn))
putRcmdr("tab6",ttkframe(tn))## called 6 because of parallel treatment with 
                             ## fractional factorial menu

tkadd(tn,tab1,text="Base Settings")   ### tabid=0
tkadd(tn,tab6,text="Export") ### tabid=5

tkconfigure(tn, takefocus=0)

nameFrame <- ttkframe(tab1)

#typeradioFrame <- ttklabelframe(tab1, text=gettextRcmdr("Augment or from scratch ?"))
#typerbVariable <- tclVar(.stored.designDopt$typerbVariable)
#augmentrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("add star portion to existing design"),variable=typerbVariable, value="augment")
#fromscratchrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("create new central composite design"),variable=typerbVariable, value="fromscratch")
#tkgrid(augmentrb, sticky="w")  ## in this case, leave default option from options
#tkgrid(fromscratchrb, sticky="w")

tkgrid(ttklabel(tab1, text=gettextRcmdr("The active data frame is the candidate data set:")), ttklabel(tab1, text=getRcmdr(".activeDataSet")), sticky="w")

baseFrame <- ttklabelframe(tab1,text=gettextRcmdr("Details for D-optimization"))

### widgets for tab1 and base frame
putRcmdr("nameVar", tclVar(.stored.designDopt$nameVar))
if (tclvalue(getRcmdr("nameVar"))=="Design.1.Dopt" & !getRcmdr(".activeDataSet")=="Design.1")
putRcmdr("nameVar", tclVar(paste(getRcmdr(".activeDataSet"), "Dopt", sep=".")))
    
nameEntry <- tkentry(nameFrame, width="20", textvariable=nameVar)
    tkbind(nameEntry, "<FocusIn>", nameenter)
    tkbind(nameEntry, "<FocusOut>", namechange)

nrunsVar <- tclVar(.stored.designDopt$nrunsVar)
nrunsEntry <- tkentry(tab1, width="8", textvariable=nrunsVar)

rhsVariable <- tclVar(.stored.designDopt$rhsVariable)
## formulaEntry is no longer used, formula is entered in frames created by modelFormulaDoE
formulaEntry <- tkentry(baseFrame, width="80", textvariable=rhsVariable)

putRcmdr("constraintVar", tclVar(.stored.designDopt$constraintVar))
constraintEntry <- tkentry(baseFrame, width="80", textvariable=constraintVar)

nrepeatVar <- tclVar(.stored.designDopt$nrepeatVar)
nrepeatEntry <- tkentry(baseFrame, width="12", textvariable=nrepeatVar)
nrepeatHint <- ttklabel(baseFrame, text="(positive integer number, large=long run time)", foreground="#888888")

randomizecbVariable <- tclVar(.stored.designDopt$cbInitials[2])
randomizecb <- ttkcheckbutton(baseFrame,text=gettextRcmdr("Randomize design"),variable=randomizecbVariable)
tkconfigure(randomizecb, takefocus=0)
seedVar <- tclVar(sample(31999,1))  ## always new
seedEntry <- tkentry(baseFrame, width="8", textvariable=seedVar)
tkconfigure(seedEntry, takefocus=0)


### Finalize tab1
helptab1Button <- buttonRcmdr(nameFrame, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab1, 
        default = "normal", borderwidth = 3)
## grid nameFrame
tkconfigure(helptab1Button, takefocus=0)
tkgrid(tklabel(nameFrame, text="Name of new design"), nameEntry, helptab1Button, sticky="w")
tkgrid(nameFrame, sticky="w", columnspan=4, pady=20)
tkgrid.configure(helptab1Button, sticky="ne", padx=15)

## grid nruns
tkgrid(nrunslab <- tklabel(tab1, text=gettextRcmdr("Numer of runs (required)")), nrunsEntry, sticky="w", pady=c(5,10))

## grid formula information
    ## define variables accessible for the formula
    Variables(colnames(get(ActiveDataSet())))
    modelFrame <- ttklabelframe(tab1, text="Specify model for optimization")
    modelFormulaDoE(hasLhs=FALSE, rhschr=tclvalue(rhsVariable))  ## creates and grids several frames 
    ### does not use rhschr; why? Is correctly available but unused
    tkgrid(modelFrame, sticky="w", columnspan="2")
    subsetBoxDoE()   ## for constraints
    tkgrid(subsetFrame, sticky="w", columnspan="2", pady=5)
    

## grid base frame
tkgrid(nRepeatlab <- tklabel(baseFrame, text=gettextRcmdr("Number of optimization repeats")), nrepeatEntry, nrepeatHint, sticky="w", pady=15)
tkgrid(randlab <- tklabel(baseFrame, text="You normally do not need to change randomization settings"),sticky="w",columnspan=3)
tkgrid(seedlab <- tklabel(baseFrame, text=gettextRcmdr("Seed for randomization")), seedEntry, 
       randomizecb, sticky="w")
tkgrid(baseFrame, sticky="nw",columnspan=3)

## tab6 for exporting
helptab6Button <- buttonRcmdr(tab6, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab6, 
        default = "normal", borderwidth = 3)

exportlabVar <- nameVar
exportlab <- ttklabel(tab6, textvariable=exportlabVar)
tkgrid(ttklabel(tab6,text="Current design to be saved:"),exportlab,helptab6Button,sticky="w") 
tkgrid.configure(exportlab, pady=15)
tkgrid.configure(helptab6Button, sticky="ne")

## radio buttons for choosing export type
etradioFrame <- ttklabelframe(tab6, text=gettextRcmdr("(How to) Export ?"))
etyperbVariable <- tclVar(.stored.designDopt$etyperbVariable)
noexprb <- tkradiobutton(etradioFrame,text=gettextRcmdr("no export"),variable=etyperbVariable,value="none")
allrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("all file types"),variable=etyperbVariable,value="all")
rdarb <- tkradiobutton(etradioFrame,text=gettextRcmdr("rda only"),variable=etyperbVariable,value="rda")
htmlrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("html and rda"),variable=etyperbVariable,value="html")
csvrb <- tkradiobutton(etradioFrame,text=gettextRcmdr("csv and rda"),variable=etyperbVariable,value="csv")
tkgrid(noexprb, sticky="w")
tkgrid(allrb, sticky="w")
tkgrid(rdarb, sticky="w")
tkgrid(htmlrb, sticky="w")
tkgrid(csvrb, sticky="w")

## radio buttons for choosing export decimal separator
decimalradioFrame <- ttklabelframe(tab6, text=gettextRcmdr("Decimal Separator ?"))
decimalrbVariable <- tclVar(.stored.designDopt$decimalrbVariable)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(tab6, text=gettextRcmdr("Storage Directory"))
putRcmdr("dirVar", tclVar(.stored.designDopt$dirVar))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")

## export file name
putRcmdr("fileVar", tclVar(.stored.designDopt$fileVar))
fileEntry <- tkentry(tab6, width="20", textvariable=fileVar)
efnamelabel <- tklabel(tab6,text=gettextRcmdr("Export file names: name below with appropriate endings (html or csv, and rda)"))
replacecbVariable <- tclVar(.stored.designDopt$cbInitials[8])
replacecb <- ttkcheckbutton(tab6,text=gettextRcmdr("Replace file(s), if exists"),variable=replacecbVariable)

## always grid details, as otherwise default file name does not work
## design name info and help button have already been gridded above
tkgrid(etradioFrame, decimalradioFrame, sticky="nw")
tkgrid(dirFrame, sticky="w", columnspan=5)
tkgrid.configure(dirFrame, pady=15)
tkgrid(efnamelabel, sticky="w", columnspan=5)
tkgrid(fileEntry, sticky="w", columnspan=5)
tkgrid(replacecb, sticky="w", columnspan=5)


## add buttons outside the notebook
buttonFrame <- tkframe(topdes2)
## die sind aber nicht dunkelgruen ...
refreshButton <- buttonRcmdr(buttonFrame, text = gettextRcmdr("Refresh form"), 
        foreground = "darkgreen", width = "12", command = onRefresh, 
        default = "normal", borderwidth = 3)
storeButton <- buttonRcmdr(buttonFrame, text = gettextRcmdr("Store form"), 
        foreground = "darkgreen", width = "12", command = onStore, 
        default = "normal", borderwidth = 3)
loadButton <- buttonRcmdr(buttonFrame, text = gettextRcmdr("Load form"), 
        foreground = "darkgreen", width = "12", command = onLoad, 
        default = "normal", borderwidth = 3)
resetButton <- buttonRcmdr(buttonFrame, text = gettextRcmdr("Reset form"), 
        foreground = "darkgreen", width = "12", command = onReset, 
        default = "normal", borderwidth = 3)
#        tkgrid(refreshButton,sticky="w")
#        tkgrid(tklabel(buttonFrame,text="  "),sticky="w")
        tkgrid(storeButton,sticky="w")
        tkgrid(loadButton,sticky="w")
        tkgrid(resetButton,sticky="w")

tkconfigure(refreshButton, takefocus=0)
tkconfigure(storeButton, takefocus=0)
tkconfigure(loadButton, takefocus=0)
tkconfigure(resetButton, takefocus=0)

## storage buttons to the right of the notebook
tkgrid(tn, buttonFrame, sticky="w", columnspan=2)

OKCancelHelp(window=topdes2, helpSubject="Menu.Dopt")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="ew")

### relations among widgets
if (exists("activestab.tn", where="RcmdrEnv")){
                tcl(tn, "select", activestab.tn)
                rm(activestab.tn, pos="RcmdrEnv")
                }

dialogSuffix(window=topdes2, rows=2, columns=2, focus=tn, bindReturn=FALSE)

}
###
# End of Menu.Dopt
###
