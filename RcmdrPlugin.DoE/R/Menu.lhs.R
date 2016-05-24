## two instances of assign replaced by justDoIt
Menu.lhs <- function(){
initializeDialogDoE(title=gettextRcmdr("Create latin hypercube sample ..."))   
     ## function initializeDialogDoE assumes topdes2 as windowname
     ## last stored top left corner for window is stored under topleft2xy
     ## onRefresh still makes window walk a little

if (exists("curindex", where="RcmdrEnv")) rm(curindex, pos="RcmdrEnv")

if (!exists(".stored.designlhs",where="RcmdrEnv")) 
           assign(".stored.designlhs", .default.designlhs,pos="RcmdrEnv")
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
     print(help("Menu.lhsTab1"))
 }
 onHelpTab2 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.FacDetails2Tab"))
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
        nrunVar=tclvalue(nrunVar),nfacVar=tclvalue(nfacVar),
        digitsVar=tclvalue(digitsVar),
        typerbVariable=tclvalue(typerbVariable),
        cbInitials = c("0", "1",
                       tclvalue(facnameAutoVariable),tclvalue(faclevelCommonVariable),
                       1,0,
                       0,tclvalue(replacecbVariable),0,
                       0
                       ),
        level1Var=tclvalue(level1Var),level2Var=tclvalue(level2Var),seedVar=tclvalue(seedVar),
        facnamlist=as.character(tclObj(facnamlist)),
        faclev1list=as.character(tclObj(faclev1list)),
        faclev2list=as.character(tclObj(faclev2list)),
        faclablist=as.character(tclObj(faclablist)),
        etyperbVariable=tclvalue(etyperbVariable),
        decimalrbVariable=tclvalue(decimalrbVariable),
        dirVar=tclvalue(dirVar), fileVar=tclvalue(fileVar))
        class(hilf) <- c("menu.designlhs","list")
        putRcmdr(".stored.designlhs",hilf)
}

onOK <- function(){
    onRefreshEnd()
    ## store entries so that users do not have to redo everything
    ## in case of stupid mistakes
    storeRcmdr()
    ## seed is not used from previously stored design
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)){
            errorCondition(window=topdes2,recall=Menu.lhs, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.lhs, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
    ###  further error messages with return to menu ?

    textfactornameslist.forcommand <- paste("factor.names=list(",paste(paste(as.character(tclObj(facnamlist)),"=c(",
                            dquote(as.character(tclObj(faclev1list))), ",",
                            dquote(as.character(tclObj(faclev2list))), ")",sep=""),
                            collapse=","),")")

    ### not yet perfect, especially NULL entries are not possible
    ### for didactic reasons distinguish between usage of default.levels and other?
    command <- paste("lhs.design(","type=",dquote(tclvalue(typerbVariable)),
                  ", nruns=",tclvalue(nrunVar),",nfactors=",tclvalue(nfacVar),
                  ",digits=",tclvalue(digitsVar), ",seed=",tclvalue(seedVar),
                  ",",textfactornameslist.forcommand,")") 

        hilf <- justDoItDoE(command)
        closeDialog(window=topdes2)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.lhs, message=gettextRcmdr(hilf))
             return()
            }
                  
        logger(paste(name, "<-", command))
        logger("## creator element of design.info will be different, when using the command line command!")
        ## change creator to contain menu settings
        hilfatt <- design.info(hilf)
        hilfatt$creator <- .stored.designlhs
        class(hilfatt$creator) <- c("menu.designlhs", "list")
        attr(hilf, "design.info") <- hilfatt
        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(name, hilf, envir=.GlobalEnv)
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
            errorCondition(window=topdes2,recall=Menu.lhs, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        }
        rm(activestab.tn, pos="RcmdrEnv")
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

listDesignlhs <- function (envir = .GlobalEnv, ...) 
{
    Vars <- ls(envir = envir, all.names = TRUE)
    Vars[which(sapply(Vars, function(.x){
               aus <- FALSE
               if ("menu.designlhs" %in% class(get(.x, envir = envir))) aus <- TRUE
               else if ("design" %in% class(get(.x, envir = envir)))
                    if ("menu.designlhs" %in% class(design.info(get(.x, envir = envir))$creator))
                       aus <- TRUE
               aus
               }))]
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
        putRcmdr(".stored.designlhs",get(lb$varlist[as.numeric(tclvalue(tcl(lb$listbox, "curselection")))+1]))
        if ("design" %in% class(getRcmdr(".stored.designlhs"))) 
            putRcmdr(".stored.designlhs", design.info(getRcmdr(".stored.designlhs"))$creator)
        tkfocus(CommanderWindow())
        tkdestroy(topdes2)
        tkdestroy(deschoose2)
        Menu.lhs()
    }
    OKCancelHelp(window=deschoose2)
    tkgrid(buttonsFrame, sticky="s")
    dialogSuffix(window=deschoose2, rows=1, columns=1, 
         focus=lb$listbox)
}

onRefreshEnd <- function(){
        nfacchange()
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
          Menu.lhs()
}

onStore <- function(){
        ## Speichernamen abfragen und hier ermöglichen (statt stored.designlhs)
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
        ## replace assign by justDoIt; assign(savename.RcmdrPlugin.DoE, getRcmdr(".stored.designlhs"), envir=.GlobalEnv)
        justDoIt(paste(savename.RcmdrPlugin.DoE, "<- getRcmdr(\".stored.designlhs\")"))
        message(gettextRcmdr("inputs have been stored"))
        }
}

onReset <- function(){
        assign(".stored.designlhs",.default.designlhs,pos="RcmdrEnv")
        tkfocus(CommanderWindow())
  tkdestroy(topdes2)
  Menu.lhs()
}

    nfacchange <- function(){
        nfacold <- length(as.character(tclObj(varlistshort)))
        nfacnew <- as.numeric(tclvalue(nfacVar))
        if (nfacold==nfacnew) return()
        if (nfacnew < nfacold){
           varlistshortt <- if (nfacnew<=50) 
                 Letters[1:nfacnew] else paste("F",1:nfacnew,sep="")
           putRcmdr("varlistshortt" , varlistshortt)
           putRcmdr("varlistshort", tclVar(getRcmdr("varlistshortt")))
           putRcmdr("facnamlist", tclVar(as.character(tclObj(facnamlist))[1:nfacnew]))
           putRcmdr("faclev1list", tclVar(as.character(tclObj(faclev1list))[1:nfacnew]))
           putRcmdr("faclev2list", tclVar(as.character(tclObj(faclev2list))[1:nfacnew]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[1:nfacnew]))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(faclev1ListBox, listvariable=faclev1list, height=min(10,nfacnew))
           tkconfigure(faclev2ListBox, listvariable=faclev2list, height=min(10,nfacnew))
           tkconfigure(faclabListBox, listvariable=faclablist, height=min(10,nfacnew))
           tkconfigure(facnameListBox, listvariable=facnamlist, height=min(10,nfacnew))
             if (selpos > nfacnew){
                tcl(fsel, "current", "0")
                factorsel()
             }
           }
        if (nfacnew > nfacold){
           varlistshortt <- if (nfacnew<=50) 
                 Letters[1:nfacnew] else paste("F",1:nfacnew,sep="")
           putRcmdr("varlistshortt" , varlistshortt)
           putRcmdr("varlistshort", tclVar(getRcmdr("varlistshortt")))
           putRcmdr("facnamlist", tclVar(c(as.character(tclObj(facnamlist)),
                               getRcmdr("varlistshortt")[(nfacold+1):nfacnew])) )
           putRcmdr("faclev1list", tclVar(c(as.character(tclObj(faclev1list)),
                               rep(tclvalue(level1Var),nfacnew-nfacold))))
           putRcmdr("faclev2list", tclVar(c(as.character(tclObj(faclev2list)),
                               rep(tclvalue(level2Var),nfacnew-nfacold))))
           putRcmdr("faclablist", tclVar(c(as.character(tclObj(faclablist)),
                               rep("",nfacnew-nfacold))))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(facnameListBox, listvariable=facnamlist, height=min(10,nfacnew))
           tkconfigure(faclev1ListBox, listvariable=faclev1list, height=min(10,nfacnew))
           tkconfigure(faclev2ListBox, listvariable=faclev2list, height=min(10,nfacnew))
           tkconfigure(faclabListBox, listvariable=faclablist, height=min(10,nfacnew))
        }
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
    factorsel<-function(){
        #### aendert die in der Textbox dargestellte Auswahl
        #### ruiniert aber leider auch wieder die korrekte Ueberschreibung der Werte
        putRcmdr("selpos", as.numeric(tclvalue(tcl(fsel, "current")))+1)
        putRcmdr("curfac", tclVar(as.character(tclObj(varlistshort))[selpos]))
        putRcmdr("curfnam", tclVar(as.character(tclObj(facnamlist))[selpos]))
        putRcmdr("curflev1", tclVar(as.character(tclObj(faclev1list))[selpos]))
        putRcmdr("curflev2", tclVar(as.character(tclObj(faclev2list))[selpos]))
        putRcmdr("curflab", tclVar(as.character(tclObj(faclablist))[selpos]))
        tkconfigure(fnam, textvariable=curfnam)
        tkconfigure(flev1, textvariable=curflev1)
        tkconfigure(flev2, textvariable=curflev2)
        tkconfigure(flab, textvariable=curflab)
    }
    fnamchange <- function(){
        ## selpos known from factorsel
        if (is.valid.name(tclvalue(curfnam))){
        hilf <- as.character(tclObj(facnamlist))
        hilf[selpos] <- tclvalue(curfnam)
        putRcmdr("facnamlist",tclVar(hilf))
        ### "facnamlist" is not automatically updated in the listbox
        ### therefore the tkconfigure
        tkconfigure(facnameListBox, listvariable=facnamlist)
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid factor name")
    }
     level1enter <- function(){
              putRcmdr("the.common.level1", tclvalue(getRcmdr("level1Var")))
        }
     level1change <- function(){
        if (identical(getRcmdr("the.common.level1"), tclvalue(getRcmdr("level1Var")))) return()
        onRefresh()
    }
     level2enter <- function(){
              putRcmdr("the.common.level2", tclvalue(getRcmdr("level2Var")))
        }
     level2change <- function(){
        if (identical(getRcmdr("the.common.level2"), tclvalue(getRcmdr("level2Var")))) return()
        onRefresh()
    }
    flev1change <- function(){
        ## selpos known from factorsel
        if (length(as.character(tclObj(curflev1)))==1){
          hilf <- as.character(tclObj(faclev1list))
          hilf[selpos] <- tclvalue(curflev1)
          putRcmdr("faclev1list",tclVar(hilf))
          tkconfigure(faclev1ListBox, listvariable=faclev1list)
        }
        else tkmessageBox(message="Empty entries or entries with blanks are not permitted, please correct!",
            icon="error", type="ok", title="Invalid factor level")
    }    
    flev2change <- function(){
        ## selpos known from factorsel
        if (length(as.character(tclObj(curflev2)))==1){
          hilf <- as.character(tclObj(faclev2list))
          hilf[selpos] <- tclvalue(curflev2)
          putRcmdr("faclev2list",tclVar(hilf))
          tkconfigure(faclev2ListBox, listvariable=faclev2list)
        }
        else tkmessageBox(message="Empty entries or entries with blanks are not permitted, please correct!",
            icon="error", type="ok", title="Invalid factor level")
    }
    flabchange <- function(){
        ## selpos known from factorsel
        ## for FocusOut event on flab
        ## still problematic, if Focus out occurs with tab
        ## as there is also a tab key event
          hilf <- as.character(tclObj(faclablist))
          hilf[selpos] <- tclvalue(curflab)
          putRcmdr("faclablist",tclVar(hilf))
          tkconfigure(faclabListBox, listvariable=faclablist)
    }
         
    tabflab <- function(){
        ## for Tab key event on flab
        ## the traversal still jumps to the first traversable control on the sheet 
        ## (rather than staying with fnam, if asked by tkfocus to do so)
        ## takefocus has so far been set to 0 for all widgets except the factor detail ones on this tab
        flabchange()  ## otherwise, not carried out!
        hilf <- as.numeric(tclvalue(tcl(fsel,"current")))+1
        if (hilf  >= as.numeric(tclvalue(nfacVar))) return()
        tcl(fsel,"current", hilf)
        factorsel()
        #tkfocus(fnam)
        #tcl(fnam, "selection", "range", 1, "end")
        #tcl("break")
    }


    swap <- function(a,b){
        hilf <- 1:as.numeric(tclvalue(nfacVar))
        hilf[b] <- a
        hilf[a] <- b
        hilf
    }

    indexchange <- function(){
        if (curindex < as.numeric(tclvalue(nfacVar)))
             putRcmdr("orderDown",swap(curindex, curindex+1))
        if (curindex > 1)
             putRcmdr("orderUp",swap(curindex, curindex-1))
        tcl(fsel, "current", curindex-1)
        factorsel()
    }

    checkIndexShort <- function(){
        putRcmdr("curindex", as.numeric(tcl(facshortListBox,"curselection"))+1)
        indexchange()
    }
    checkIndexNam <- function(){
        putRcmdr("curindex", as.numeric(tcl(facnameListBox,"curselection"))+1)
        indexchange()
    }
    checkIndexLev1 <- function(){
        putRcmdr("curindex", as.numeric(tcl(faclev1ListBox,"curselection"))+1)
        indexchange()
    }
    checkIndexLev2 <- function(){
        putRcmdr("curindex", as.numeric(tcl(faclev2ListBox,"curselection"))+1)
        indexchange()
    }
    checkIndexLab <- function(){
        putRcmdr("curindex", as.numeric(tcl(faclabListBox,"curselection"))+1)
        indexchange()
    }


    onUp <- function(){
        if (!exists("curindex")) return()
        if (length(curindex)==0) return()
        if (curindex=="1" | is.null(curindex)) return()
        else {
           putRcmdr("facnamlist", tclVar(as.character(tclObj(facnamlist))[orderUp]))
           putRcmdr("faclev1list", tclVar(as.character(tclObj(faclev1list))[orderUp]))
           putRcmdr("faclev2list", tclVar(as.character(tclObj(faclev2list))[orderUp]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderUp]))
           tkconfigure(faclev1ListBox, listvariable=faclev1list)
           tkconfigure(faclev2ListBox, listvariable=faclev2list)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex-1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(faclev1ListBox,"selection","set",curindex-1)
           tcl(faclev2ListBox,"selection","set",curindex-1)
           tcl(faclabListBox,"selection","set",curindex-1)
           tcl(facnameListBox,"selection","set",curindex-1)
           }
    }

    onDown <- function(){
        if (!exists("curindex")) return()
        if (length(curindex)==0) return()
        if (curindex==as.numeric(tclvalue(nfacVar)) | is.null(curindex)) return()
        else {
           putRcmdr("facnamlist", tclVar(as.character(tclObj(facnamlist))[orderDown]))
           putRcmdr("faclev1list", tclVar(as.character(tclObj(faclev1list))[orderDown]))
           putRcmdr("faclev2list", tclVar(as.character(tclObj(faclev2list))[orderDown]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderDown]))
           tkconfigure(faclev1ListBox, listvariable=faclev1list)
           tkconfigure(faclev2ListBox, listvariable=faclev2list)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex+1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(faclev1ListBox,"selection","set",curindex-1)
           tcl(faclev2ListBox,"selection","set",curindex-1)
           tcl(faclabListBox,"selection","set",curindex-1)
           tcl(facnameListBox,"selection","set",curindex-1)
           }
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
        }
 }

######## end define functions                          


##### define userform
#tn <- ttknotebook(top,height=100, width=500)

putRcmdr("tn",ttknotebook(topdes2))
#tn <- ttknotebook(topdes2)

putRcmdr("tab1",ttkframe(tn))
putRcmdr("tab2",ttkframe(tn))
putRcmdr("tab6",ttkframe(tn))## called 6 because of parallel treatment with 
                             ## fractional factorial menu

tkadd(tn,tab1,text="Base Settings")   ### tabid=0
tkadd(tn,tab2,text="Factor Details")  ### tabid=1
tkadd(tn,tab6,text="Export") ### tabid=5

tkconfigure(tn, takefocus=0)

nameFrame <- ttkframe(tab1)
typeradioFrame <- ttklabelframe(tab1, text=gettextRcmdr("Type of lhs ?"))
baseFrame <- ttklabelframe(tab1,text=gettextRcmdr("Size and randomization"))

### widgets for tab1 and base frame
putRcmdr("nameVar", tclVar(.stored.designlhs$nameVar))
nameEntry <- tkentry(nameFrame, width="20", textvariable=nameVar)
    tkbind(nameEntry, "<FocusIn>", nameenter)
    tkbind(nameEntry, "<FocusOut>", namechange)

typerbVariable <- tclVar(.stored.designlhs$typerbVariable)
optimumrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("optimum"),variable=typerbVariable, value="optimum")
geneticrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("genetic"),variable=typerbVariable, value="genetic")
improvedrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("improved"),variable=typerbVariable, value="improved")
maximinrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("maximin"),variable=typerbVariable, value="maximin")
randomrb <- tkradiobutton(typeradioFrame,text=gettextRcmdr("random"),variable=typerbVariable, value="random")
tkgrid(optimumrb, sticky="w")  
tkgrid(geneticrb, sticky="w")
tkgrid(improvedrb, sticky="w")
tkgrid(maximinrb, sticky="w")
tkgrid(randomrb, sticky="w")

nrunVar <- tclVar(.stored.designlhs$nrunVar)
nrunEntry <- tkentry(baseFrame, width="8", textvariable=nrunVar)
nrunHint <- ttklabel(baseFrame, text="(positive integer number)", foreground="#888888")
nfacVar <- tclVar(.stored.designlhs$nfacVar)
nfacEntry <- tkentry(baseFrame, width="8", textvariable=nfacVar)
nfacHint <- ttklabel(baseFrame, text="(positive integer number)", foreground="#888888")
tkbind(nfacEntry,"<FocusOut>",nfacchange)
digitsVar <- tclVar(.stored.designlhs$digitsVar)
digitsEntry <- tkentry(baseFrame, width="8", textvariable=digitsVar)
digitsHint <- ttklabel(baseFrame, text="(integer number)", foreground="#888888")

## radio buttons for choosing design type
seedVar <- tclVar(sample(31999,1))  ## always new
seedEntry <- tkentry(baseFrame, width="8", textvariable=seedVar)
seedHint <- tklabel(baseFrame, text="(You normally do not need to change this)", foreground="#888888")
tkconfigure(seedEntry, takefocus=0)

## preparations for bottom frame
bottomFrame <- tkframe(topdes2)

## grid base frame
tkgrid(nrunlab <- tklabel(baseFrame, text=gettextRcmdr("Number of runs")), nrunEntry, nrunHint, sticky="w")
## omitted nfaccb, on form, nfactors must always be specified
tkgrid(nfaclab <- tklabel(baseFrame, text=gettextRcmdr("Number of factors")), nfacEntry, nfacHint, sticky="w")
tkgrid(digitslab <- tklabel(baseFrame, text=gettextRcmdr("Number of decimal places")), 
    digitsEntry, digitsHint, sticky="w", pady="15")
tkgrid(seedlab <- tklabel(baseFrame, text=gettextRcmdr("Seed for randomization")), seedEntry, sticky="w")
tkgrid(seedHint, sticky="w")

helptab1Button <- buttonRcmdr(nameFrame, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab1, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab1Button, takefocus=0)

### Finalize tab1
tkgrid(tklabel(nameFrame, text="Name of new design"), nameEntry, helptab1Button, sticky="w")
tkgrid(nameFrame, sticky="w", columnspan=4)
tkgrid.configure(nameFrame, pady=40)
tkgrid.configure(helptab1Button, sticky="ne")
tkgrid(typeradioFrame, baseFrame, sticky="nw",columnspan=3)
tkgrid.configure(baseFrame, padx=10)

## Factor Details Tab
## factor details frame
### facnameAutoVariable (not needed any more) and faclevelCommonVariable

## default levels frame
deflevFrame <- ttklabelframe(tab2,text="Default levels (lower and upper scale ends)")
facnameAutoVariable <- tclVar(.stored.designlhs$cbInitials[3])
faclevelCommonVariable <- tclVar(.stored.designlhs$cbInitials[4])
faclevelCommonButton <- ttkcheckbutton(deflevFrame,text=gettextRcmdr("Common factor levels"),
    variable=faclevelCommonVariable,command=onRefresh)
tkconfigure(faclevelCommonButton,takefocus=0)
putRcmdr("level1Var", tclVar(.stored.designlhs$level1Var))
    level1Entry <- ttkentry(deflevFrame, width="20", textvariable=level1Var)
    tkconfigure(level1Entry,takefocus=0)
    tkbind(level1Entry, "<FocusIn>", level1enter)
    tkbind(level1Entry, "<FocusOut>", level1change)
tkconfigure(level1Entry,takefocus=0)
putRcmdr("level2Var", tclVar(.stored.designlhs$level2Var))
    level2Entry <- tkentry(deflevFrame, width="20", textvariable=level2Var)
    tkconfigure(level2Entry,takefocus=0)
    tkbind(level2Entry, "<FocusIn>", level2enter)
    tkbind(level2Entry, "<FocusOut>", level2change)
    tkgrid(faclevelCommonButton,sticky="w",columnspan=3)
faclevCommonLab<-tklabel(deflevFrame,text=gettextRcmdr("CAUTION: Checking this box overwrites all custom factor levels."))
if (!as.logical(as.numeric(tclvalue(faclevelCommonVariable)))){ 
    tkgrid(faclevCommonLab,sticky="w", columnspan=3)
    tkgrid.configure(faclevCommonLab,pady=10)
    }
tkgrid(tklabel(deflevFrame, text=gettextRcmdr("Lower scale end")),tklabel(deflevFrame,text="  ",width=2),tklabel(deflevFrame, text=gettextRcmdr("Upper scale end")),sticky="e")
tkgrid(level1Entry, tklabel(deflevFrame,text="  ",width=2),level2Entry, sticky="e")

## factor details
## values as vectors
facnamlistt <- .stored.designlhs$facnamlist
if (as.logical(as.numeric(tclvalue(faclevelCommonVariable)))) {
    faclev1listt <- rep(tclvalue(level1Var),tclvalue(nfacVar)) 
    faclev2listt <- rep(tclvalue(level2Var),tclvalue(nfacVar)) 
    } else{
    faclev1listt <- .stored.designlhs$faclev1list
    faclev2listt <- .stored.designlhs$faclev2list
    }
faclablistt <- .stored.designlhs$faclablist
varlistshortt <- if (as.numeric(tclvalue(nfacVar))<=50) 
                 Letters[1:tclvalue(nfacVar)] else paste("F",1:tclvalue(nfacVar),sep="")

    enterlistFrame <- ttkframe(tab2)
    listFrame <- ttklabelframe(enterlistFrame, text="Factor Details")
    putRcmdr("selpos", 1)
    putRcmdr("curfac", tclVar(varlistshortt[1]))
    putRcmdr("curfnam", tclVar(facnamlistt[1]))
    putRcmdr("curflev1", tclVar(faclev1listt[1]))
    putRcmdr("curflev2", tclVar(faclev2listt[1]))
    putRcmdr("curflab", tclVar(faclablistt[1]))
    
        ## fsel must select the right factor
    ## this should be highlighted in factor lists
    ##    and all related entries shown for changing in text boxes fnam etc.
    enterFrame <- ttklabelframe(enterlistFrame, text=gettextRcmdr("Modify factor details for selected factor"))
    fsel <- ttkcombobox(enterFrame, textvariable=curfac, width=5, values=varlistshortt, state="readonly")
    tkbind(fsel, "<<ComboboxSelected>>", factorsel)
    #fnam <- ttkentry(listFrame, textvariable=curfnam, width=20,validate="focusout", validatecommand=fnamchange)
    fnam <- ttkentry(enterFrame, textvariable=curfnam, width=15)
    tkbind(fnam, "<FocusOut>", fnamchange)
    flev1 <- ttkentry(enterFrame, textvariable=curflev1, width=15)
    tkbind(flev1, "<FocusOut>", flev1change)
    if (as.logical(as.numeric(tclvalue(faclevelCommonVariable)))){ 
        tkconfigure(flev1,state="disabled")
        }
    flev2 <- ttkentry(enterFrame, textvariable=curflev2, width=15)
    tkbind(flev2, "<FocusOut>", flev2change)
    if (as.logical(as.numeric(tclvalue(faclevelCommonVariable)))){
        tkconfigure(flev2,state="disabled")
        }
    flab <- ttkentry(enterFrame, textvariable=curflab, width=20)
    tkbind(flab, "<FocusOut>", flabchange)
    tkbind(flab, "<Key-Tab>", tabflab)
    tkgrid(tklabel(enterFrame,text=gettextRcmdr("Select"),width=6),
           tklabel(enterFrame,text=gettextRcmdr("Factor name"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("First level"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("Second level"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("Comment or label \n(for html export only)"), width=20),
           sticky="w")
    tkgrid(fsel,fnam, flev1, flev2, flab, sticky="w")
    
    putRcmdr("facnamlist", tclVar(facnamlistt))
    putRcmdr("varlistshort", tclVar(varlistshortt))
    putRcmdr("faclev1list", tclVar(faclev1listt))
    putRcmdr("faclev2list", tclVar(faclev2listt))
    putRcmdr("faclablist", tclVar(faclablistt))

    facshortListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=varlistshort,
        width = 6, background="#EBEBDC")
    tkbind(facshortListBox, "<<TraverseIn>>",function() tkfocus(fsel))

    facnameListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=facnamlist,
        width = 15, background="#EBEBDC")
    faclev1ListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclev1list,
        width = 15, background="#EBEBDC")
    faclev2ListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclev2list,
        width = 15, background="#EBEBDC")
    faclabListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclablist,
        width = 20, background="#EBEBDC")

    ## determine current index and reordering for onUp and onDown
    tkbind(facshortListBox, "<<ListboxSelect>>", checkIndexShort)
    tkbind(facnameListBox, "<<ListboxSelect>>", checkIndexNam)
    tkbind(faclev1ListBox, "<<ListboxSelect>>", checkIndexLev1)
    tkbind(faclev2ListBox, "<<ListboxSelect>>", checkIndexLev2)
    tkbind(faclabListBox, "<<ListboxSelect>>", checkIndexLab)


    ### funktioniert, ist aber noch nicht schön
    scrollbar <- ttkscrollbar(listFrame, command = function(...) {
            tkyview(facshortListBox, ...)
            tkyview(facnameListBox, ...)
            tkyview(faclev1ListBox, ...)
            tkyview(faclev2ListBox, ...)
            tkyview(faclabListBox, ...)
            })

#    tkgrid(tklabel(enterlistFrame,text="  ", width=5),enterFrame, sticky="w")
    tkgrid(enterFrame, sticky="w", columnspan=5)
    tkgrid.configure(enterFrame, pady=10)
    ## Hoch-/Runterschieben von Einträgen ermöglichen

    downupFrame <- ttkframe(listFrame)
    moveDownButton <- buttonRcmdr(downupFrame, text = gettextRcmdr("Move Down"), 
        foreground = "darkgreen", command = onDown, 
        default = "normal", borderwidth = 3, width=12)
    moveUpButton <- buttonRcmdr(downupFrame, text = gettextRcmdr("Move Up"), 
        foreground = "darkgreen", command = onUp, 
        default = "normal", borderwidth = 3, width=12)
    tkgrid(moveDownButton, sticky="w")
    tkgrid(moveUpButton, sticky="w")

    tkgrid(scrollbar, facshortListBox, facnameListBox, faclev1ListBox, faclev2ListBox, faclabListBox, downupFrame, sticky = "nw")
    tkgrid.configure(scrollbar, sticky = "wns")
    tkgrid.configure(facnameListBox, sticky = "new")

helptab2Button <- buttonRcmdr(tab2, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab2, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab2Button, takefocus=0)

## finalize tab2 Factor details
    tkgrid(helptab2Button, sticky="e")
    tkgrid(deflevFrame, sticky="nw")
    tkgrid.configure(deflevFrame, pady=10)
    tkgrid(listFrame, columnspan=6,sticky="w")
    tkgrid(enterlistFrame, columnspan=6,sticky="w")

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
etyperbVariable <- tclVar(.stored.designlhs$etyperbVariable)
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
decimalrbVariable <- tclVar(.stored.designlhs$decimalrbVariable)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(tab6, text=gettextRcmdr("Storage Directory"))
putRcmdr("dirVar", tclVar(.stored.designlhs$dirVar))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")

## export file name
putRcmdr("fileVar", tclVar(.stored.designlhs$fileVar))
fileEntry <- tkentry(tab6, width="20", textvariable=fileVar)
efnamelabel <- tklabel(tab6,text=gettextRcmdr("Export file names: name below with appropriate endings (html or csv, and rda)"))
replacecbVariable <- tclVar(.stored.designlhs$cbInitials[8])
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

OKCancelHelp(window=topdes2, helpSubject="Menu.lhs")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, bottomFrame, sticky="ew")


### relations among widgets
if (exists("activestab.tn", where="RcmdrEnv")){
                tcl(tn, "select", activestab.tn)
                rm(activestab.tn, pos="RcmdrEnv")
                }

dialogSuffix(window=topdes2, rows=2, columns=2, focus=tn, bindReturn=FALSE)

}
###
# End of Menu.lhs
###
