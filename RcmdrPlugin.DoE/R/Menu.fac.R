## two instances of assign replaced by justDoIt

Menu.fac <- function(){
initializeDialogDoE(title=gettextRcmdr("Create full factorial design ..."))   
     ## function initializeDialogDoE assumes topdes2 as windowname
     ## last stored top left corner for window is stored under topleft2xy
     ## onRefresh still makes window walk a little
     
     ## Lesen der level klappt noch nicht, jedenfalls Fehler dass nlevels und factor.names nicht zusammenpassen
     
     ## menus must be edited

if (exists("curindex", where="RcmdrEnv")) rm(curindex, pos="RcmdrEnv")

if (!exists(".stored.designfac",where="RcmdrEnv")) 
           assign(".stored.designfac", .default.designfac,pos="RcmdrEnv")
           ## nameVar, nrunVar, nfacVar, nrepVar, nblockVar, 
           ## cbInitials containing repeat.onlyVariable, randomizeVariable, 
           ##                       facnamesAutoVariable, faclevelsCommonVariable, 
           ##                       nrunVar, estcbVariable
           ##                       specialcbVariable, replacecbVariable, MaxC2cbVariable
           ##                       res3cbVariable
           ## nlevVar, level2Var, seedVar, specialrbVariable, hardVar, genVar, 
           ## catlgVar, designVar, designrbVariable, destyperbVariable
           ## resVar, qualcritrbVariable, facnamlist,nlevlist,faclevlist, faclablist
           ## etyperbVariable, decimalrbVariable, dirVar, fileVar
## support change o implementing nblockVar
           if (is.null(getRcmdr(".stored.designfac"))) 
               putRcmdr(".stored.designfac", 
               c(list(nblockVar="1"), getRcmdr(".stored.designfac")))

## MaxC2cbVariable is free again (no. 9 of cbInitials)

## define called functions
 infoClose <- function(){
     putRcmdr("infotxt",tclVar(""))
 }
 
 onHelpTab1 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.facTab1"))
 }
 onHelpTab2 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.FacDetailsGenTab"))
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
        nrunVar=tclvalue(nrunVar),nfacVar=tclvalue(nfacVar),nrepVar=tclvalue(nrepVar), 
        nblockVar=tclvalue(nblockVar), 
        cbInitials = c(tclvalue(repeat.onlyVariable), tclvalue(randomizeVariable),
                       0,0,
                       1,0,
                       0,tclvalue(replacecbVariable),0,
                       0
                       ),
        seedVar=tclvalue(seedVar),
        facnamlist=as.character(tclObj(facnamlist)),
        nlevlist=as.character(tclObj(nlevlist)),
        faclevlist=as.character(tclObj(faclevlist)),
        faclablist=as.character(tclObj(faclablist)),
        etyperbVariable=tclvalue(etyperbVariable),
        decimalrbVariable=tclvalue(decimalrbVariable),
        dirVar=tclvalue(dirVar), fileVar=tclvalue(fileVar))
        class(hilf) <- c("menu.designfac","list")
        putRcmdr(".stored.designfac",hilf)
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
            errorCondition(window=topdes2,recall=Menu.fac, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Object")))){
              errorCondition(window=topdes2,recall=Menu.fac, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
        if (any(getRcmdr(".stored.designfac")$nlevlist=="") )
          {
              errorCondition(window=topdes2,recall=Menu.fac, 
              gettextRcmdr("Factor details are not completely specified."))
              return()
          }
        if (any(getRcmdr(".stored.designfac")$faclevlist=="") )
          {
              errorCondition(window=topdes2,recall=Menu.fac, 
              gettextRcmdr("Factor details are not completely specified."))
              return()
          }
    ###  further error messages with return to menu ?

    textfactornameslist.forcommand <- paste("factor.names=list(",paste(paste(as.character(tclObj(facnamlist)),"=c(",
                            sapply(strsplit(as.character(tclObj(faclevlist))," "),function(obj) paste(dquote(obj),collapse=",")), 
                            ")",sep=""),collapse=","),")")

    ### not yet perfect, especially NULL entries are not possible
    command <- paste("fac.design(nfactors=",tclvalue(nfacVar),",replications=",
                  tclvalue(nrepVar),",repeat.only=",as.logical(as.numeric(tclvalue(repeat.onlyVariable))),
                  ",blocks=",tclvalue(nblockVar),
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",nlevels=c(", paste(as.character(tclObj(nlevlist)),collapse=","),
                  "),",textfactornameslist.forcommand,")") 

                  
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.fac, message=gettextRcmdr(hilf))
             return()
            }
        logger(paste(name, "<-", command))
        logger("## creator element of design.info will be different, when using the command line command!")
        ## change creator to contain menu settings
        hilfatt <- design.info(hilf)
        hilfatt$creator <- .stored.designfac
        class(hilfatt$creator) <- c("menu.designfac", "list")
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
                stop("file ", paste(filename, "rda", "."), " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
          if (tclvalue(etyperbVariable)=="html" & tolower(paste(filename, "html", sep = ".")) %in% lf) 
                stop("file ", paste(filename, "html", "."), " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
          if (tclvalue(etyperbVariable)=="csv" & tolower(paste(filename, "csv", sep = ".")) %in% lf) 
                stop("file ", paste(filename, "csv", "."), " exists and must not be replaced. Change filename on Export tab or allow replacing of files.")
         }
        if (tclvalue(decimalrbVariable)=="default") command <- paste("export.design(",name,
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(path),", file=",dquote(filename),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),")",sep="")
        else command <- paste("export.design(",name, 
               ", type=",dquote(tclvalue(etyperbVariable)),",path=",dquote(path),", file=",dquote(filename),", replace=",
               as.logical(as.numeric(tclvalue(replacecbVariable))),", OutDec=", dquote(tclvalue(decimalrbVariable)),")",sep="")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            errorCondition(window=topdes2,recall=Menu.fac, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        }
        rm(activestab.tn, pos="RcmdrEnv")
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

listDesignfac <- function (envir = .GlobalEnv, ...) 
{
    Vars <- ls(envir = envir, all.names = TRUE)
    Vars[which(sapply(Vars, function(.x){
               aus <- FALSE
               if ("menu.designfac" %in% class(get(.x, envir = envir))) aus <- TRUE
               else if ("design" %in% class(get(.x, envir = envir)))
                    if ("menu.designfac" %in% class(design.info(get(.x, envir = envir))$creator))
                       aus <- TRUE
               aus
               }))]
}


onLoad <- function(){
    ## seems to work now, needs to be tested!
        hilf <- listDesignfac()
        if (length(hilf)==0) {
            tkmessageBox(message=gettextRcmdr("There are no stored design inputs in this session."),icon="error", type="ok", title="no stored design inputs")
            return()
            }
    putRcmdr("deschoosefac",tktoplevel())
    tkwm.title(deschoosefac, gettextRcmdr("Choose stored design form"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(deschoosefac, position)
    putRcmdr("lb", variableListBox(deschoosefac, variableList=hilf, title="Choose stored design form"))
        tkgrid(lb$frame)
    onOK <- function() {
        putRcmdr(".stored.designfac", 
             get(lb$varlist[as.numeric(tclvalue(tcl(lb$listbox, "curselection")))+1]))
        if ("design" %in% class(getRcmdr(".stored.designfac"))) 
            putRcmdr(".stored.designfac", design.info(getRcmdr(".stored.designfac"))$creator)
        if (is.null(getRcmdr(".stored.designfac")$nblockVar)) 
               putRcmdr(".stored.designfac", 
               c(list(nblockVar="1"), getRcmdr(".stored.designfac")))
        tkfocus(CommanderWindow())
        tkdestroy(topdes2)
        tkdestroy(deschoosefac)
        Menu.fac()
    }
    OKCancelHelp(window=deschoosefac)
    tkgrid(buttonsFrame, sticky="s")
    dialogSuffix(window=deschoosefac, rows=1, columns=1, 
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
          Menu.fac()
}

onStore <- function(){
        ## Speichernamen abfragen und hier ermöglichen (statt stored.designfac)
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
        ## replace assign by justDoIt; assign(savename.RcmdrPlugin.DoE, getRcmdr(".stored.designfac"), envir=.GlobalEnv)
        justDoIt(paste(savename.RcmdrPlugin.DoE, "<- getRcmdr(\".stored.designfac\")"))
        message(gettextRcmdr("inputs have been stored"))
        }
}

onReset <- function(){
        assign(".stored.designfac",.default.designfac,pos="RcmdrEnv")
        tkfocus(CommanderWindow())
  tkdestroy(topdes2)
  Menu.fac()
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
           putRcmdr("nlevlist", tclVar(as.character(tclObj(nlevlist))[1:nfacnew]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[1:nfacnew]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[1:nfacnew]))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(nlevListBox, listvariable=nlevlist, height=min(10,nfacnew))
           tkconfigure(faclevListBox, listvariable=faclevlist, height=min(10,nfacnew))
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
           putRcmdr("nlevlist", tclVar(c(as.character(tclObj(nlevlist)),
                               rep("",nfacnew-nfacold))))
           putRcmdr("faclevlist", tclVar(c(as.character(tclObj(faclevlist)),
                               rep("",nfacnew-nfacold))))
           putRcmdr("faclablist", tclVar(c(as.character(tclObj(faclablist)),
                               rep("",nfacnew-nfacold))))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(facnameListBox, listvariable=facnamlist, height=min(10,nfacnew))
           tkconfigure(nlevListBox, listvariable=nlevlist, height=min(10,nfacnew))
           tkconfigure(faclevListBox, listvariable=faclevlist, height=min(10,nfacnew))
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
        putRcmdr("curnlev", tclVar(as.character(tclObj(nlevlist))[selpos]))
        putRcmdr("curflev", tclVar(as.character(tclObj(faclevlist))[selpos]))
        putRcmdr("curflab", tclVar(as.character(tclObj(faclablist))[selpos]))
        tkconfigure(fnam, textvariable=curfnam)
        tkconfigure(nlev, textvariable=curnlev)
        tkconfigure(flev, textvariable=curflev)
        tkconfigure(flab, textvariable=curflab)
    }
    fnamchange <- function(){
        ## selpos known from factorsel
        if (is.valid.name(tclvalue(curfnam))){
        hilf <- as.character(tclObj(facnamlist))
        hilf[selpos] <- tclvalue(curfnam)
        ### "facnamlist" is not automatically updated in the listbox
        ### therefore the tkconfigure
        putRcmdr("facnamlist",tclVar(hilf))
        tkconfigure(facnameListBox, listvariable=facnamlist)
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid factor name")
    }
    nlevchange <- function(){
        ## selpos known from factorsel
        if (length(as.character(tclObj(curnlev)))==1){
          hilf <- as.character(tclObj(nlevlist))
          hilf[selpos] <- tclvalue(curnlev)
          putRcmdr("nlevlist",tclVar(hilf))
          tkconfigure(nlevListBox, listvariable=nlevlist)
        }
        else tkmessageBox(message="Empty entries or entries with blanks are not permitted, please correct!",
            icon="error", type="ok", title="Invalid number of factor levels")
        if (is.na(as.numeric(tclvalue(curnlev))) | !floor(as.numeric(tclvalue(curnlev)))==as.numeric(tclvalue(curnlev))) 
            tkmessageBox(message="an integer number was expected, please correct!",
            icon="error", type="ok", title="Invalid number of levels")
        nrunVar <- tclVar(prod(as.numeric(hilf)))
        tkconfigure(nrunShow, textvariable=nrunVar)
        if (length(as.character(tclObj(curflev)))==0){
            flevchange()
            tkconfigure(flev, textvariable=curflev)}
    }    
    flevchange <- function(){
        ## selpos known from factorsel
        nlev <- as.numeric(tclvalue(curnlev))
        nlevspec <- length(as.character(tclObj(curflev)))
        if (nlevspec==nlev){
          hilf <- as.character(tclObj(faclevlist))
          hilf[selpos] <- tclvalue(curflev)
          ### updating hilf does work
          ### but "varlist" is not automatically updated in the listbox
          ### therefore the tkconfigure
          putRcmdr("faclevlist",tclVar(hilf))
          tkconfigure(faclevListBox, listvariable=faclevlist)
        }
        else{ 
            if (nlevspec==0){ 
              putRcmdr("curflev", tclVar(paste(1:nlev, sep=" ")))
              hilf <- as.character(tclObj(faclevlist))
              hilf[selpos] <- tclvalue(curflev)
              putRcmdr("faclevlist",tclVar(hilf))
              tkconfigure(faclevListBox, listvariable=faclevlist)
            }
            else tkmessageBox(message=gettextRcmdr(nlev," levels entries needed, \nyou specified ", nlevspec, " levels, please correct!"),
                icon="error", type="ok", title="Invalid number of level entries")
        }
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
    checkIndexNlev <- function(){
        putRcmdr("curindex", as.numeric(tcl(nlevListBox,"curselection"))+1)
        indexchange()
    }
    checkIndexFlev <- function(){
        putRcmdr("curindex", as.numeric(tcl(faclevListBox,"curselection"))+1)
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
           putRcmdr("nlevlist", tclVar(as.character(tclObj(nlevlist))[orderUp]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[orderUp]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderUp]))
           tkconfigure(nlevListBox, listvariable=nlevlist)
           tkconfigure(faclevListBox, listvariable=faclevlist)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex-1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(nlevListBox,"selection","set",curindex-1)
           tcl(faclevListBox,"selection","set",curindex-1)
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
           putRcmdr("nlevlist", tclVar(as.character(tclObj(nlevlist))[orderDown]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[orderDown]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderDown]))
           tkconfigure(nlevListBox, listvariable=nlevlist)
           tkconfigure(faclevListBox, listvariable=faclevlist)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex+1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(nlevListBox,"selection","set",curindex-1)
           tcl(faclevListBox,"selection","set",curindex-1)
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
baseFrame <- ttklabelframe(tab1,text=gettextRcmdr("Size and randomization"))

### widgets for tab1 and base frame
putRcmdr("nameVar", tclVar(.stored.designfac$nameVar))
nameEntry <- tkentry(nameFrame, width="20", textvariable=nameVar)
    tkbind(nameEntry, "<FocusIn>", nameenter)
    tkbind(nameEntry, "<FocusOut>", namechange)

nrunVar <- tclVar(.stored.designfac$nrunVar)
nrunShow <- ttklabel(baseFrame, width="8", textvariable=nrunVar)
nrunHint <- ttklabel(baseFrame, text="(product of all numbers of factor levels)", foreground="#888888")
nfacVar <- tclVar(.stored.designfac$nfacVar)
nfacEntry <- tkentry(baseFrame, width="8", textvariable=nfacVar)
tkbind(nfacEntry,"<FocusOut>",nfacchange)
nrepVar <- tclVar(.stored.designfac$nrepVar)
nrepEntry <- tkentry(baseFrame, width="8", textvariable=nrepVar)
nblockVar <- tclVar(.stored.designfac$nblockVar)
randomizeVariable <-  tclVar(.stored.designfac$cbInitials[2])
nblockEntry <- tkentry(baseFrame, width="8", textvariable=nblockVar)
randomizecb <- ttkcheckbutton(baseFrame,text=gettextRcmdr("Randomization"),variable=randomizeVariable)
tkconfigure(randomizecb, takefocus=0)
seedVar <- tclVar(sample(31999,1))  ## always new
seedEntry <- tkentry(baseFrame, width="8", textvariable=seedVar)
tkconfigure(seedEntry, takefocus=0)
repeat.onlyVariable <- tclVar(.stored.designfac$cbInitials[1])
repeat.onlycb <- ttkcheckbutton(baseFrame,text=gettextRcmdr("Repeat only"),variable=repeat.onlyVariable)
tkconfigure(repeat.onlycb, takefocus=0)

## grid base frame
tkgrid(nrunlab <- tklabel(baseFrame, text=gettextRcmdr("Number of runs")), nrunShow, nrunHint, sticky="w")
## omitted nfaccb, on form, nfactors must always be specified
tkgrid(nfaclab <- tklabel(baseFrame, text=gettextRcmdr("Number of factors")), nfacEntry, sticky="w")
tkgrid(nblocklab <- tklabel(baseFrame, text=gettextRcmdr("Number of blocks")), nblockEntry, sticky="w")
tkgrid(nreplab <- tklabel(baseFrame, text=gettextRcmdr("Replications")), nrepEntry, repeat.onlycb, sticky="w")
tkgrid.configure(nreplab, pady=15)
tkgrid(randlab <- tklabel(baseFrame, text="You normally do not need to change randomization settings"),sticky="w",columnspan=3)
tkgrid(seedlab <- tklabel(baseFrame, text=gettextRcmdr("Seed for randomization")), seedEntry, 
       randomizecb, sticky="w")

helptab1Button <- buttonRcmdr(nameFrame, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab1, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab1Button, takefocus=0)

### Finalize tab1
tkgrid(tklabel(nameFrame, text="Name of new design"), nameEntry, helptab1Button, sticky="w")
tkgrid(nameFrame, sticky="w", columnspan=4)
tkgrid.configure(nameFrame, pady=40)
tkgrid.configure(helptab1Button, sticky="ne")
tkgrid(baseFrame, sticky="nw",columnspan=3)


## Factor Details Tab
## factor details frame
### facnameAutoVariable (not needed any more) and faclevelCommonVariable

## factor details
## values as vectors
    facnamlistt <- .stored.designfac$facnamlist
    nlevlistt <- .stored.designfac$nlevlist
    faclevlistt <- .stored.designfac$faclevlist
    faclablistt <- .stored.designfac$faclablist
    varlistshortt <- if (as.numeric(tclvalue(nfacVar))<=50) 
                 Letters[1:tclvalue(nfacVar)] else paste("F",1:tclvalue(nfacVar),sep="")

    enterlistFrame <- ttkframe(tab2)
    listFrame <- ttklabelframe(enterlistFrame, text="Factor Details")
    putRcmdr("selpos", 1)
    putRcmdr("curfac", tclVar(varlistshortt[1]))
    putRcmdr("curfnam", tclVar(facnamlistt[1]))
    putRcmdr("curnlev", tclVar(nlevlistt[1]))
    putRcmdr("curflev", tclVar(faclevlistt[1]))
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
    nlev <- ttkentry(enterFrame, textvariable=curnlev, width=6)
    tkbind(nlev, "<FocusOut>", nlevchange)
    flev <- ttkentry(enterFrame, textvariable=curflev, width=20)
    tkbind(flev, "<FocusOut>", flevchange)
    flab <- ttkentry(enterFrame, textvariable=curflab, width=20)
    tkbind(flab, "<FocusOut>", flabchange)
    tkbind(flab, "<Key-Tab>", tabflab)
    tkgrid(tklabel(enterFrame,text=gettextRcmdr("Select"),width=6),
           tklabel(enterFrame,text=gettextRcmdr("Factor name"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("No. of \nlevels"), width=6),
           tklabel(enterFrame,text=gettextRcmdr("Levels \n(separate with blanks)"), width=20),
           tklabel(enterFrame,text=gettextRcmdr("Comment or label \n(for html export only)"), width=20),
           sticky="w")
    tkgrid(fsel,fnam, nlev, flev, flab, sticky="w")
    
    putRcmdr("facnamlist", tclVar(facnamlistt))
    putRcmdr("varlistshort", tclVar(varlistshortt))
    putRcmdr("nlevlist", tclVar(nlevlistt))
    putRcmdr("faclevlist", tclVar(faclevlistt))
    putRcmdr("faclablist", tclVar(faclablistt))

    facshortListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=varlistshort,
        width = 6, background="#EBEBDC")
    tkbind(facshortListBox, "<<TraverseIn>>",function() tkfocus(fsel))

    facnameListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=facnamlist,
        width = 15, background="#EBEBDC")
    nlevListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=nlevlist,
        width = 6, background="#EBEBDC")
    faclevListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclevlist,
        width = 20, background="#EBEBDC")
    faclabListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclablist,
        width = 20, background="#EBEBDC")

    ## determine current index and reordering for onUp and onDown
    tkbind(facshortListBox, "<<ListboxSelect>>", checkIndexShort)
    tkbind(facnameListBox, "<<ListboxSelect>>", checkIndexNam)
    tkbind(nlevListBox, "<<ListboxSelect>>", checkIndexNlev)
    tkbind(faclevListBox, "<<ListboxSelect>>", checkIndexFlev)
    tkbind(faclabListBox, "<<ListboxSelect>>", checkIndexLab)


    ### funktioniert, ist aber noch nicht schön
    scrollbar <- ttkscrollbar(listFrame, command = function(...) {
            tkyview(facshortListBox, ...)
            tkyview(facnameListBox, ...)
            tkyview(nlevListBox, ...)
            tkyview(faclevListBox, ...)
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

    tkgrid(scrollbar, facshortListBox, facnameListBox, nlevListBox, faclevListBox, faclabListBox, downupFrame, sticky = "nw")
    tkgrid.configure(scrollbar, sticky = "wns")
    tkgrid.configure(facnameListBox, sticky = "new")

helptab2Button <- buttonRcmdr(tab2, text = gettextRcmdr("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab2, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab2Button, takefocus=0)

## finalize tab2 Factor details
    tkgrid(helptab2Button, sticky="e")
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
etyperbVariable <- tclVar(.stored.designfac$etyperbVariable)
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
decimalrbVariable <- tclVar(.stored.designfac$decimalrbVariable)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(tab6, text=gettextRcmdr("Storage Directory"))
putRcmdr("dirVar", tclVar(.stored.designfac$dirVar))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")

## export file name
putRcmdr("fileVar", tclVar(.stored.designfac$fileVar))
fileEntry <- tkentry(tab6, width="20", textvariable=fileVar)
efnamelabel <- tklabel(tab6,text=gettextRcmdr("Export file names: name below with appropriate endings (html or csv, and rda)"))
replacecbVariable <- tclVar(.stored.designfac$cbInitials[8])
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

OKCancelHelp(window=topdes2, helpSubject="Menu.fac")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, sticky="ew")


### relations among widgets
if (!as.logical(as.numeric(tclvalue(randomizeVariable)))){
tkconfigure(seedEntry, state="disabled")
tkconfigure(seedlab, state="disabled")
}else {
tkconfigure(seedEntry, state="normal")
tkconfigure(seedlab, state="normal")
}
if (exists("activestab.tn", where="RcmdrEnv")){
                tcl(tn, "select", activestab.tn)
                rm(activestab.tn, pos="RcmdrEnv")
                }

dialogSuffix(window=topdes2, rows=2, columns=2, focus=tn, bindReturn=FALSE)

}
###
# End of Menu.fac
###
