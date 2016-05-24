## two instances of assign replaced by justDoItDoE with modified command and justDoIt("design.info(name) <- hilfatt")
## one instance of assign replaced by justDoIt

Menu.FrF2level <- function(){

initializeDialogDoE(title=gettext("Create regular 2-level design ..."))   
     ## function initializeDialogDoE assumes topdes2 as windowname
     ## last stored top left corner for window is stored under topleft2xy
     ## onRefresh still makes window walk a little

if (exists("curindex", where="RcmdrEnv")) rm(curindex, pos="RcmdrEnv")

if (!exists(".stored.design2FrF", where="RcmdrEnv")) 
           assign(".stored.design2FrF", .default.design2,pos="RcmdrEnv")
           ## nameVar, nrunVar, nfacVar, nrepVar
           ## cbInitials containing repeat.onlyVariable, randomizeVariable, 
           ##                       aliasblock2fiVariable, faclevelsCommonVariable, 
           ##                       nrunEntryVariable, estcbVariable
           ##                       specialcbVariable, replacecbVariable, MaxC2cbVariable
           ##                       res3cbVariable
           ## level1Var, level2Var, seedVar, specialrbVariable, hardVar, genVar, 
           ## catlgVar, designVar, designrbVariable, destyperbVariable
           ## resVar, qualcritrbVariable, facnamlist,faclev1list,faclev2list, faclablist
           ## estrbVariable, maxtimeVar, est2fislist,
           ## etyperbVariable, decimalrbVariable, dirVar, fileVar

## MaxC2cbVariable is free again (no. 9 of cbInitials)

## define called functions
 infoClose <- function(){
     putRcmdr("infotxt",tclVar(""))
 }
 
 onHelpTab1 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.FrF2levelTab1"))
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
 
 onHelpTabEstimable <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows")
            tkgrab.release(topdes2)
     print(help("Menu.FrF2levelTabEstimable"))
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
        nblockVar=tclvalue(nblockVar), ncenterVar=tclvalue(ncenterVar),
        cbInitials = c(tclvalue(repeat.onlyVariable), tclvalue(randomizeVariable),
                       tclvalue(aliasblock2fiVariable),tclvalue(faclevelCommonVariable),
                       tclvalue(nrunEntryVariable),0,
                       tclvalue(specialcbVariable),tclvalue(replacecbVariable),0,
                       tclvalue(res3cbVariable)
                       ),
        level1Var=tclvalue(level1Var),level2Var=tclvalue(level2Var),seedVar=tclvalue(seedVar),
        specialrbVariable=tclvalue(specialrbVariable),hardVar=tclvalue(hardVar),
        designrbVariable=tclvalue(designrbVariable),
        genVar=tclvalue(genVar),catlgVar=tclvalue(catlgVar),designVar=tclvalue(designVar),
        resVar=tclvalue(resVar),
        qualcritrbVariable=tclvalue(qualcritrbVariable),
        comprclassVar=tclvalue(comprclassVar),
        facnamlist=as.character(tclObj(facnamlist)),
        faclev1list=as.character(tclObj(faclev1list)),
        faclev2list=as.character(tclObj(faclev2list)),
        faclablist=as.character(tclObj(faclablist)),
        estrbVariable=tclvalue(estrbVariable),
        comprrbVariable=tclvalue(comprrbVariable),
        maxtimeVar=tclvalue(maxtimeVar),est2fislist=est2fislist,
        etyperbVariable=tclvalue(etyperbVariable),
        decimalrbVariable=tclvalue(decimalrbVariable),
        dirVar=tclvalue(dirVar), fileVar=tclvalue(fileVar))
    class(hilf) <- c("menu.design2FrF","list")
    putRcmdr(".stored.design2FrF", hilf)
}

onOK <- function(){
    onRefreshEnd()
    ## store entries so that users do not have to redo everything
    ## in case of stupid mistakes
    ## seed is not used from previously stored design
    storeRcmdr()
     closeDialog(window=topdes2)
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(window=topdes2,recall=Menu.FrF2level, 
                    message=paste('"', name, '" ', gettext("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(name, gettext("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.FrF2level, 
              message=gettext("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }

    ###  further error messages with return to menu ?

    textfactornameslist.forcommand <- paste("factor.names=list(",paste(paste(as.character(tclObj(facnamlist)),"=c(",
                            dquote(as.character(tclObj(faclev1list))), ",",
                            dquote(as.character(tclObj(faclev2list))), ")",sep=""),
                            collapse=","),")")

    ### not yet perfect, especially NULL entries are not possible
    ### also, not very didactical, as default settings are unnecessarily included
    ### 
    
    ## do always
        MaxC2 <- FALSE
        if (tclvalue(qualcritrbVariable)=="MaxC2") MaxC2 <- TRUE
        nrun.forcommand <- tclvalue(nrunVar)
        resolution <- "NULL"
        alias.block.2fis <- "FALSE"
        if (as.logical(as.numeric(as.character(tclvalue(aliasblock2fiVariable)))))
            alias.block.2fis <- TRUE
    if (!as.logical(as.numeric(as.character(tclvalue(nrunEntryVariable))))){
          nrun.forcommand <- "NULL"
          resolution <- 3
          if (tclvalue(resVar)=="IV") resolution <- 4
          if (tclvalue(resVar)=="V+") resolution <- 5
        }
    
    if (!as.logical(as.numeric(as.character(tclvalue(specialcbVariable))))){
            ### separate statements because of repeat.only replications in blocked designs
            if (!(as.logical(as.numeric(tclvalue(repeat.onlyVariable))) & as.numeric(tclvalue(nblockVar))>0)) 
            command <- paste("FrF2(nruns=",nrun.forcommand,",nfactors=",tclvalue(nfacVar),", blocks=", tclvalue(nblockVar),
                  ", alias.block.2fis =", alias.block.2fis, 
                  ", ncenter=", tclvalue(ncenterVar), ", MaxC2 =", MaxC2, ", resolution =", resolution,
                  ",replications=",tclvalue(nrepVar),",repeat.only=",as.logical(as.numeric(tclvalue(repeat.onlyVariable))),
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",",textfactornameslist.forcommand,")")
            else 
            command <- paste("FrF2(nruns=",nrun.forcommand,",nfactors=",tclvalue(nfacVar),", blocks=", tclvalue(nblockVar),
                  ", alias.block.2fis =", alias.block.2fis, 
                  ", ncenter=", tclvalue(ncenterVar), ", MaxC2 =", MaxC2, ", resolution =", resolution,
                  ",wbreps=",tclvalue(nrepVar),",repeat.only=TRUE, randomize=",
                  as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",",textfactornameslist.forcommand,")")
                  }
    else{
               estimable <- "NULL"
          ## only do if special cases present
          if (!tclvalue(estrbVariable)=="none" & length(est2fislist)>0){
                ## only if selected and specified!!
                ### estimable interactions
                ## compromise plans
                if (tclvalue(comprrbVariable) == "compr") {
                    command <- paste("compromise(",tclvalue(nfacVar),", c(",
                       paste(dquote(which(Letters  %in% notest2fislist)),collapse=","),
                       "), ",substr(tclvalue(comprclassVar),1,1),")")
                    hilf <- justDoItDoE(command)
                    if (class(hilf)[1]=="try-error") {
                          Message(paste(gettext("Offending command:"), "\n", command), type="error")
                          errorCondition(window=topdes2,recall=Menu.FrF2level, message=gettext(hilf))
                          return()
                    }
                    ## replace assign by justDoIt; assign("calc.estim", hilf, envir=.GlobalEnv)
                    putRcmdr("hilf", hilf)
                    justDoIt("calc.estim <- getRcmdr(\"hilf\")")
                    rm("hilf", pos="RcmdrEnv")

                    #compromise(as.numeric(tclvalue(nfacVar)), which(Letters  %in% notest2fislist), 
                    #        as.numeric(substr(tclvalue(comprclassVar),1,1))))
                    logger(paste("calc.estim <-", command))
                    estimable <- "calc.estim$requirement"
                    if (tclvalue(estrbVariable)=="distinct") estimable <- paste(estimable, ", perms=calc.estim$perms.full")
                }
                else estimable <- paste("c(",paste(dquote(est2fislist),collapse=","),")")
                
                if (!(tclvalue(hardVar)=="0" & tclvalue(designrbVariable)=="default" & tclvalue(nblockVar)=="1"))
                    tk_messageBox(message=gettext("estimable has taken precedence, not all other requests have been granted!"),type="ok")
               clear <- "TRUE"
               if (tclvalue(estrbVariable)=="distinct") clear <- "FALSE"
               res3 <- as.logical(as.numeric(as.character(tclvalue(res3cbVariable))))
            command <- paste("FrF2(nruns=",nrun.forcommand,",nfactors=",tclvalue(nfacVar), ", MaxC2 =", MaxC2, 
                  ",replications=",tclvalue(nrepVar),",repeat.only=",as.logical(as.numeric(tclvalue(repeat.onlyVariable))),
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",",textfactornameslist.forcommand,
                  ", estimable=", estimable, ", clear =", clear, ", res3 =", res3, ", max.time =", tclvalue(maxtimeVar), 
                  ", select.catlg =", tclvalue(catlgVar),       ")")
          }
          else {
            hard.forcommand <- tclvalue(hardVar)
            if (hard.forcommand == "0") hard.forcommand <- "NULL"
            generators <- "NULL"
            design <- "NULL"
            if (tclvalue(designrbVariable)=="gen") generators <- paste("c(",paste(dquote(unlist(strsplit(tclvalue(genVar),","))),collapse=","),")")
            if (tclvalue(designrbVariable)=="design") design <- dquote(tclvalue(designVar))
            
            command <- paste("FrF2(nruns =",nrun.forcommand,",nfactors =",tclvalue(nfacVar),
                  ", blocks =", tclvalue(nblockVar), ", alias.block.2fis =", alias.block.2fis, 
                  ", ncenter =", tclvalue(ncenterVar), 
                  ", hard=",hard.forcommand, ", generators =", generators, 
                  ", design =", design, 
                  ",replications=",tclvalue(nrepVar),",repeat.only=",as.logical(as.numeric(tclvalue(repeat.onlyVariable))),
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",",textfactornameslist.forcommand, ", select.catlg =", tclvalue(catlgVar), ")")
          }
        }       ## end of special           
        hilf <- justDoItDoE(command)
        if (tclvalue(estrbVariable)=="distinct" & length(est2fislist)>0){
                 diagcommand <- paste("print(",dQuote(paste(gettext("Design search in progress: you allowed up to"), 
                        tclvalue(maxtimeVar), gettext("seconds"))), ")")
                 doItAndPrint(diagcommand, log=FALSE)
                 }
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettext("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.FrF2level, message=gettext(hilf))
#            if (tclvalue(estrbVariable)=="distinct" & length(est2fislist)>0){
#                 diagcommand <- paste("print(",dQuote(paste(gettext("No design found in"), 
#                       tclvalue(maxtimeVar), gettext("seconds"))), ")")
#                 doItAndPrint(diagcommand, log=FALSE)
#                 diagcommand <- paste("print(",dQuote(gettext("Experts may try to speed up the search using command line programming (?estimable.2fis).")), ")")
#                 doItAndPrint(diagcommand, log=FALSE)
#                 }
             return()
            }
        logger(paste(name, "<-", command))
        logger("## creator element of design.info will be different, when using the command line command!")
        ## change creator to contain menu settings
        hilfatt <- design.info(hilf)
        hilfatt$creator <- .stored.design2FrF
        class(hilfatt$creator) <- c("menu.design2FrF", "list")
        
        attr(hilf, "design.info") <- hilfatt
        ## replace assign by justDoIt; assign(name, hilf, envir=.GlobalEnv)
        putRcmdr("hilf", hilf)
        justDoIt(paste(name, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
        activeDataSet(name)
        ## remove calc.estim
        if (as.logical(as.numeric(as.character(tclvalue(specialcbVariable)))) & tclvalue(comprrbVariable)=="compr"){ 
            rm(calc.estim, envir=.GlobalEnv)
            logger("rm(calc.estim)")
        }
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
            errorCondition(window=topdes2,recall=Menu.FrF2level, message=gettext(hilf))
             return()
            }
        logger(command)
        }
        rm(activestab.tn, pos="RcmdrEnv")
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

listDesign2 <- function (envir = .GlobalEnv, ...) 
{
    Vars <- ls(envir = envir, all.names = TRUE)
    Vars[which(sapply(Vars, function(.x){
               aus <- FALSE
               if ("menu.design2FrF" %in% class(get(.x, envir = envir))) aus <- TRUE
               else if ("design" %in% class(get(.x, envir = envir)))
                    if ("menu.design2FrF" %in% class(design.info(get(.x, envir = envir))$creator))
                       aus <- TRUE
               aus
               }))]
}


onLoad <- function(){
    ## seems to work now, needs to be tested!
        hilf <- listDesign2()
        if (length(hilf)==0) {
            tkmessageBox(message=gettext("There are no stored design inputs in this session."),icon="error", type="ok", title="no stored design inputs")
            return()
            }
    putRcmdr("deschoose2",tktoplevel())
    tkwm.title(deschoose2, gettext("Choose stored design form"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(deschoose2, position)
    putRcmdr("lb", variableListBox(deschoose2, variableList=hilf, title="Choose stored design form"))
        tkgrid(lb$frame)
    onOK <- function() {
        putRcmdr(".stored.design2FrF",get(lb$varlist[as.numeric(tclvalue(tcl(lb$listbox, "curselection")))+1]))
        if ("design" %in% class(getRcmdr(".stored.design2FrF"))) 
            putRcmdr(".stored.design2FrF", design.info(getRcmdr(".stored.design2FrF"))$creator)
        tkfocus(CommanderWindow())
        tkdestroy(topdes2)
        tkdestroy(deschoose2)
        Menu.FrF2level()
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
        onRefreshEnd()
        ## letzte Position enthaelt tab index (beginnend bei 1)
          tkfocus(CommanderWindow())
          tkdestroy(topdes2)
          Menu.FrF2level()
}

onestrb <- function(){
        onestrb.worefresh()
        onRefresh()
}

onestrb.worefresh <- function(){
        if (!tclvalue(estrbVariable)=="none"){
             tkconfigure(selectButton, state="normal")
             tkconfigure(deselectButton, state="normal")
             tkconfigure(comprestrb, state="normal")
             tkconfigure(manualestrb, state="normal")
             tkconfigure(comprclassEntry, state="normal")
             }
        else {
             tkconfigure(selectButton, state="disabled")
             tkconfigure(deselectButton, state="disabled")
             tkconfigure(comprestrb, state="disabled")
             tkconfigure(manualestrb, state="disabled")
             tkconfigure(comprclassEntry, state="disabled")
        }
}


oncomprestrb <- function(){
        oncomprestrb.worefresh
        onRefresh()
}
oncomprestrb.worefresh <- function(){
        if (tclvalue(comprrbVariable)=="compr"){
             tkconfigure(comprclassEntry, state="normal")
             }
        else {
             tkconfigure(comprclassEntry, state="disabled")
        }
}

onSpecialcb <- function(){
       if (tclvalue(specialcbVariable)=="0"){
            if (tabpos() %in% c(2,3,4)) tcl(tn,"select",0)
           putRcmdr("estrbVariable", tclVar("none"))
           putRcmdr("specialrbVariable", tclVar("none"))
           putRcmdr("hardVar", tclVar("0"))
           putRcmdr("genVar", tclVar("NULL"))
           putRcmdr("catlgVar", tclVar("catlg"))
           putRcmdr("designVar", tclVar("NULL"))
           putRcmdr("estimable2fis", "")
           putRcmdr("estrbVariable", tclVar("none"))
           }
       tkconfigure(noestrb, variable=estrbVariable)
       tkconfigure(clearrb, variable=estrbVariable)
       tkconfigure(distinctrb, variable=estrbVariable)
       tkconfigure(defaultrb, variable=designrbVariable)
       tkconfigure(genrb, variable=designrbVariable)
       tkconfigure(catlgrb, variable=designrbVariable)
       tkconfigure(designrb, variable=designrbVariable)
       tkconfigure(genEntry, textvariable=genVar)
       tkconfigure(catlgEntry, textvariable=catlgVar)
       tkconfigure(designEntry, textvariable=designVar)
      # tkconfigure(nonerb, variable=specialrbVariable)
      # tkconfigure(hardrb, variable=specialrbVariable)
      # tkconfigure(debarrb, variable=specialrbVariable)
       onRefresh()
}

onStore <- function(){
        ## Speichernamen abfragen und hier ermöglichen (statt stored.design2)
        textentry() ## creates text string stored in savename.RcmdrPlugin.DoE
        if (!is.null(savename.RcmdrPlugin.DoE)){
        if (!is.valid.name(savename.RcmdrPlugin.DoE)) {
            textcorrect(gettext("This is not a valid name. Please correct:"))
            return()
          }
        if (is.element(savename.RcmdrPlugin.DoE, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(savename.RcmdrPlugin.DoE, gettext("Object"))))
            {
              textcorrect(gettext("Please enter a new name:"))
              return()
             }
          }
        storeRcmdr()
        ## replace assign by justDoIt; assign(savename.RcmdrPlugin.DoE, getRcmdr(".stored.design2FrF"), envir=.GlobalEnv)
        justDoIt(paste(savename.RcmdrPlugin.DoE, "<- getRcmdr(\".stored.design2FrF\")"))
        message(gettext("inputs have been stored"))
        }
}

onReset <- function(){
        assign(".stored.design2FrF",.default.design2,pos="RcmdrEnv")
        tkfocus(CommanderWindow())
  tkdestroy(topdes2)
  Menu.FrF2level()
}

    nfacchange <- function(){
        nfacold <- length(as.character(tclObj(varlistshort)))
        nfacnew <- as.numeric(tclvalue(nfacVar))
        if (nfacold==nfacnew) return()
        if (as.logical(as.numeric(as.character(tclvalue(nrunEntryVariable)))))
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors in"), tclvalue(nrunVar), gettext("runs\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))
        else 
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))
        tkconfigure(infoButton1, textvariable=infoknopftext)
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
           if (tclvalue(specialcbVariable)=="1") onRefresh()
           }
        if (nfacnew > nfacold){
           varlistshortt <- if (nfacnew<=50) 
                 Letters[1:nfacnew] else paste("F",1:nfacnew,sep="")
           putRcmdr("varlistshortt" , varlistshortt)
           putRcmdr("varlistshort", tclVar(getRcmdr("varlistshortt")))
           putRcmdr("facnamlist",tclVar(c(as.character(tclObj(facnamlist)),
                               getRcmdr("varlistshortt")[(nfacold+1):nfacnew])))
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
           if (tclvalue(specialcbVariable)=="1") onRefresh()
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
#          putRcmdr("exportlabVar", tclVar(paste("Current design to be saved:", tclvalue(nameVar),"\n   ")))  ## otherwise, variables would be directly tied
          tkconfigure(fileEntry, textvariable=getRcmdr("fileVar"))
#          tkconfigure(exportlab, textvariable=getRcmdr("exportlabVar"))
          }
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid design name")
    }

     catlgnamechange <- function(){
        if (is.valid.name(tclvalue(catlgVar))){
          if (!exists(tclvalue(catlgVar))) 
              tkmessageBox(message="catalogue does not exist",icon="error", type="ok", title="Specified catalogue does not exist")
          else{
              if (!"catlg" %in% class(get(tclvalue(catlgVar)))) 
                  tkmessageBox(message="is not a design catalogue of class catlg",icon="error", type="ok", title="invalid catalogue specified")
              else onRefresh()
#              else { pick <- nfac(get(tclvalue(catlgVar)))==as.numeric(tclvalue(nfacVar))
#                    if (tclvalue(resVar)=="IV") pick <- pick & res(get(tclvalue(catlgVar))) > 3 else
#                    if (tclvalue(resVar)=="V+") pick <- pick & res(get(tclvalue(catlgVar))) > 4
#                    if (tclvalue(nrunEntryVariable)=="1") pick <- pick & nruns(get(tclvalue(catlgVar))) == as.numeric(tclvalue(nrunVar))
#                    if (any(pick)) putRcmdr("catlgliste", names(get(tclvalue(catlgVar)))[pick])
#                    else putRcmdr("catlgliste", "")
#                    tkconfigure(designEntry, values=catlgliste)
#              }
          }
        }
        else tkmessageBox(message="invalid name!",icon="error", type="ok", title="Invalid design name")
    }

     nrunnow <- function(){
           putRcmdr("nrunOld", tclvalue(nrunVar))
        }
     nrunchange <- function(){
        if (!tclvalue(nrunVar)==nrunOld){
            if(!2^round(log2(as.numeric(tclvalue(nrunVar))))==as.numeric(tclvalue(nrunVar))){
                 tk_messageBox(caption="invalid run number",message = gettext("invalid run number, must be power of 2"), type = "ok")
                 return()}
        if (as.logical(as.numeric(as.character(tclvalue(nrunEntryVariable)))))
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors in"), tclvalue(nrunVar), gettext("runs\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))
        else 
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))
        tkconfigure(infoButton1, textvariable=infoknopftext)
            #onRefresh()
          }
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
          ### updating hilf does work
          ### but "varlist" is not automatically updated in the listbox
          ### therefore the tkconfigure
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


onInfo <- function(){
   if (tclvalue(qualcritrbVariable)=="MA") onMAinfo()
   else onMaxC2info()
}

onInfo1 <- function(){
     ## all designs for nfacVar factors in nrunVar runs
     resmin <- 3
     if (tclvalue(resVar)=="IV") resmin <- 4
     if (tclvalue(resVar)=="V+") resmin <- 5
     if (!as.logical(as.numeric(as.character(tclvalue(nrunEntryVariable)))))
          nrun <- "all"
          else nrun <- tclvalue(nrunVar)
     command <- paste("print(",tclvalue(catlgVar),", nruns =",dquote(nrun),
        ", nfactors =", tclvalue(nfacVar), 
        ", show.alias = TRUE, MaxC2 =", as.logical(tclvalue(qualcritrbVariable)=="MaxC2"), 
        ", res.min =", resmin, ")")
     hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            if (!as.numeric(tclvalue(nfacVar))<=round(log2(as.numeric(tclvalue(nrunVar)))))
              logger("# There is no such design in the catalogue. See button Available designs.")
            if (as.numeric(tclvalue(nfacVar))==round(log2(as.numeric(tclvalue(nrunVar))))) 
              logger("# This is a full factorial design.")
            if (as.numeric(tclvalue(nfacVar))<round(log2(as.numeric(tclvalue(nrunVar))))) 
              logger("# This is a replicated full factorial design.")
            return()
            }
      ## logger(command)  ## braucht man nur für justDoIt
      tkgrab.release(topdes2)
      doItAndPrint(command)
}

onInfo2 <- function(){
     ## generators for current designs
     command="dc()"
        ## logger(command)  ## braucht man nur für justDoIt
        doItAndPrint(command)
}

onInfo3 <- function(){
     ## alias pattern for current designs
     command="dc()"
        ## logger(command)  ## braucht man nur für justDoIt
        doItAndPrint(command)
}

curdes <- function(){
     ## find current design
}

onMAinfo <- function(){
    putRcmdr("info.window",tktoplevel())
    tkwm.title(info.window, gettext("Requested information"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(info.window, position)

    tcinfo <- ttklabel(info.window,text="You must CLOSE this window for continuing work with the design dialogue.")
    tkgrid(tcinfo)
    tcl("image",  "create", "photo", "resolution.image", 
         file=system.file( "images", "resolutionimage.gif", package = "RcmdrPlugin.DoE" ) )
    tc <- tklabel(info.window, image="resolution.image")
    tkgrid(tc)
    #catlginfo <- ttklabel(info.window,text="Additional regular designs from this menu:")
    #catlginfo1 <- ttklabel(info.window,text="32 runs (up to 31 factors, res. III), 64 runs (up to 32 factors res. IV, 33 to 63 factors res. III),")
    #catlginfo2 <- ttklabel(info.window,text="128 runs (up to 64 factors, res. IV, 65 to 127 factors, res. III)")
    catlginfo3 <- ttklabel(info.window,text="Additional irregular resolution V designs from general orthogonal arrays menu:")
    catlginfo4 <- ttklabel(info.window,text="128 runs (up to 15 factors), 256 runs (up to 19 factors), 2048 runs (up to 63 factors)")

    #tkgrid(catlginfo,sticky="w")
    #tkgrid(catlginfo1)
    #tkgrid(catlginfo2)
    tkgrid(catlginfo3, sticky="w")
    tkgrid(catlginfo4)
    
      dialogSuffix(window=info.window, rows=1, columns=1, 
         focus=tc,
         onOK=function() "")
 }

onMaxC2info <- function(){
    if (tclvalue(resVar)=="III") tcl("image",  "create", "photo", "MaxC2.image", 
         file=system.file( "images", "MaxC2res3image.gif", package = "RcmdrPlugin.DoE" )  )
    else if (tclvalue(resVar)=="IV") tcl("image",  "create", "photo", "MaxC2.image", 
         file=system.file( "images", "MaxC2res4image.gif", package = "RcmdrPlugin.DoE" ))
    else {tkmessageBox(message="All 2fis are clear for designs of resolution V or more.")
          return()}
    putRcmdr("info.window",tktoplevel())
    tkwm.title(info.window, gettext("Requested information"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(info.window, position)


    tcinfo <- ttklabel(info.window,text="You must CLOSE this window for continuing work with the design dialogue.")
    tkgrid(tcinfo)
    if (tclvalue(resVar)=="III") tkgrid(ttklabel(info.window,text="NOTE: red: resolution III, yellow: resolution IV, green: resolution V+"))
    if (tclvalue(resVar)=="IV") tkgrid(ttklabel(info.window,text="NOTE: yellow: resolution IV, green: resolution V+"))
    tc <- tklabel(info.window, image="MaxC2.image")
    tkgrid(tc)
    if (tclvalue(resVar)=="III") catlginfo <- ttklabel(info.window,text=gettext("Cell entries for resolution III and IV designs:"))
    else catlginfo <- ttklabel(info.window,text=gettext("Cell entries for resolution IV designs:"))
    catlginfo2 <- ttklabel(info.window,text=gettext("maximum number of clear 2fis and maximum number of factors with all 2fis clear."))
    tkgrid(catlginfo)
    tkgrid(catlginfo2)

      dialogSuffix(window=info.window, rows=1, columns=1, 
         focus=tc,
         onOK=function() "")
 }
 
 onChangeDir <- function(){
     putRcmdr("direct",tclvalue(tkchooseDirectory()))
     if (!direct=="") {
        putRcmdr("dirVar", tclVar(direct))
        tkconfigure(dirEntry, textvariable = dirVar)
        }
 }
 
 onSelect <- function(){
     if (tclvalue(tcl(notest2fis$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- notest2fislist[as.numeric(unlist(strsplit(tclvalue(tcl(notest2fis$listbox, "curselection")), " ")))+1]
     putRcmdr("notest2fislist", setdiff(notest2fislist,add))
     putRcmdr("est2fislist", intaclistt[intaclistt %in% c(est2fislist,add)])
     add <- NULL
     tcl(notest2fis$listbox, "selection", "clear", "0", "999")
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist
 }
 onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(est2fis$listbox, "curselection"))=="") return()
     add <- est2fislist[as.numeric(unlist(strsplit(tclvalue(tcl(est2fis$listbox, "curselection")), " ")))+1]
     putRcmdr("est2fislist", setdiff(est2fislist,add))
     putRcmdr("notest2fislist", intaclistt[intaclistt %in% c(notest2fislist,add)])
     add <- NULL
     tcl(est2fis$listbox, "selection", "clear", "0", "999")
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist
 }
######## end define functions                          


##### define userform
#tn <- ttknotebook(top,height=100, width=500)

putRcmdr("tn",ttknotebook(topdes2))
#tn <- ttknotebook(topdes2)

putRcmdr("tab1",ttkframe(tn))
putRcmdr("tab2",ttkframe(tn))
putRcmdr("tab3",ttkframe(tn))
putRcmdr("tab4",ttkframe(tn))
putRcmdr("tab5",ttkframe(tn))
putRcmdr("tab6",ttkframe(tn))

tkadd(tn,tab1,text="Base Settings")   ### tabid=0
tkadd(tn,tab2,text="Factor Details")  ### tabid=1
tkadd(tn,tab3,text="Estimable Model") ### tabid=2
tkadd(tn,tab4,text="Block & Split-Plot") ### tabid=3
tkadd(tn,tab5,text="Special") ### tabid=4
tkadd(tn,tab6,text="Export") ### tabid=5

tkconfigure(tn, takefocus=0)

## define canvas for showing possibilities
### image file location must be changed!!!
### widgets for tab1 and base frame
putRcmdr("nameVar", tclVar(.stored.design2FrF$nameVar))
nameEntry <- tkentry(tab1, width="20", textvariable=nameVar)
    tkbind(nameEntry, "<FocusIn>", nameenter)
    tkbind(nameEntry, "<FocusOut>", namechange)
### wird pbcbVariable noch gebraucht ?
##pbcbVariable <- tclVar(.stored.design2FrF$cbInitials[8])
#pbcb <- ttkcheckbutton(desinfoFrame,text=gettext("Screening design (Plackett-Burman)"),variable=pbcbVariable)

baseFrame <- ttklabelframe(tab1,text=gettext("Size and randomization"))

nrunVar <- tclVar(.stored.design2FrF$nrunVar)
nrunEntry <- tkentry(baseFrame, width="8", textvariable=nrunVar)
    tkbind(nrunEntry,"<FocusIn>",nrunnow)
    tkbind(nrunEntry,"<FocusOut>",nrunchange)
nfacVar <- tclVar(.stored.design2FrF$nfacVar)
nfacEntry <- tkentry(baseFrame, width="8", textvariable=nfacVar)
tkbind(nfacEntry,"<FocusOut>",nfacchange)
ncenterVar <- tclVar(.stored.design2FrF$ncenterVar)
ncenterEntry <- tkentry(baseFrame, width="8", textvariable=ncenterVar)
nblockVar <- tclVar(.stored.design2FrF$nblockVar)
nblockEntry <- tkentry(baseFrame, width="8", textvariable=nblockVar)
aliasblock2fiVariable <- tclVar(.stored.design2FrF$cbInitials[3])
aliasblock2fi <- ttkcheckbutton(baseFrame,text=gettext("blocks may be \naliased with 2fis"),variable=aliasblock2fiVariable)
tkconfigure(aliasblock2fi, takefocus=0)
nrepVar <- tclVar(.stored.design2FrF$nrepVar)
nrepEntry <- tkentry(baseFrame, width="8", textvariable=nrepVar)
randomizeVariable <-  tclVar(.stored.design2FrF$cbInitials[2])
randomizecb <- ttkcheckbutton(baseFrame,text=gettext("Randomization"),variable=randomizeVariable)
tkconfigure(randomizecb, takefocus=0)
seedVar <- tclVar(sample(31999,1))  ## always new
seedEntry <- tkentry(baseFrame, width="8", textvariable=seedVar)
tkconfigure(seedEntry, takefocus=0)
nrunEntryVariable <- tclVar(.stored.design2FrF$cbInitials[5])
nruncb <- ttkcheckbutton(baseFrame,text=gettext("Specify nruns"),variable=nrunEntryVariable, command=onRefresh)
tkconfigure(nruncb, takefocus=0)
repeat.onlyVariable <- tclVar(.stored.design2FrF$cbInitials[1])
repeat.onlycb <- ttkcheckbutton(baseFrame,text=gettext("Repeat only"),variable=repeat.onlyVariable)
tkconfigure(repeat.onlycb, takefocus=0)

despropFrame <- ttklabelframe(tab1,text="Design properties")
descritFrame <- ttkframe(despropFrame)
desinfoFrame <- ttkframe(despropFrame)

resVar <- tclVar(.stored.design2FrF$resVar)
#resEntry <- tkentry(despropframe, textvariable=resVar)
    resEntry <- ttkcombobox(descritFrame, textvariable=resVar, width=5, values=c("III","IV","V+"), state="readonly")
    tkbind(resEntry, "<<ComboboxSelected>>", onRefresh)
 
qualcritrbVariable <- tclVar(.stored.design2FrF$qualcritrbVariable)
MArb <- tkradiobutton(descritFrame,text=gettext("MA (Maximum resolution and minimum aberration)"),variable=qualcritrbVariable,value="MA")
MaxC2rb <- tkradiobutton(descritFrame,text=gettext("MaxC2 (Maximum number of clear 2fis)"),variable=qualcritrbVariable,value="MaxC2")

## grid descritFrame
tkgrid(resEntry, tklabel(descritFrame,text=gettext("Minimum resolution\nNOTE: affects design generation\nfor MaxC2 choice\nOR unspecified number of runs only"),justify="left"))
#tkgrid(tklabel(descritFrame,text=gettext("NOTE: affects design generation for MaxC2 choice")), columnspan=2, sticky="w")
#tkgrid(tklabel(descritFrame,text=gettext("OR unspecified number of runs only")), columnspan=2, sticky="w")

tkgrid(MArb, columnspan=2, sticky="w")
tkgrid(MaxC2rb, columnspan=2, sticky="w")
tkgrid(tklabel(descritFrame,text=gettext("  ")), columnspan=2)

        if (as.logical(as.numeric(as.character(tclvalue(nrunEntryVariable)))))
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors in"), tclvalue(nrunVar), gettext("runs\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))
        else 
        putRcmdr("infoknopftext", tclVar(paste(gettext("Show best 10 designs for"), tclvalue(nfacVar), 
          gettext("factors\n     The menu remains open, \n     fetch it back after looking at designs"),sep=" ")))

infoButton <- buttonRcmdr(desinfoFrame, text = gettext("Show available designs"), 
        foreground = "darkgreen", command = onInfo, 
        default = "normal", borderwidth = 3)
infoButton1 <- buttonRcmdr(desinfoFrame, textvariable = infoknopftext, 
        foreground = "darkgreen", command = onInfo1, 
        default = "normal", borderwidth = 3)
infoButton2 <- buttonRcmdr(desinfoFrame, text = gettext("Show generators for current design"), 
        foreground = "darkgreen", command = onInfo2, 
        default = "normal", borderwidth = 3)
infoButton3 <- buttonRcmdr(desinfoFrame, text = gettext("Show alias pattern for current design"), 
        foreground = "darkgreen", command = onInfo3, 
        default = "normal", borderwidth = 3)

## grid desinfoFrame
tkgrid(infoButton,sticky="we", columnspan=2)
tkgrid(infoButton1,sticky="we", columnspan=2)
#tkgrid(infoButton2,sticky="we", columnspan=2)
#tkgrid(infoButton3,sticky="we", columnspan=2)

## grid despropFrame
tkgrid(descritFrame, desinfoFrame, sticky="w", columnspan=3)

## grid design frame
## (for expert choices)
    designFrame <- ttklabelframe(tab1,text=gettext("Expert choices"))
  ### widgets for design frame

  catlab <- tklabel(designFrame, text=gettext("Catalogue of designs:"))
  putRcmdr("catlgVar", tclVar(.stored.design2FrF$catlgVar))
  ## set to valid default, if invalid entry
  if (!exists(tclvalue(catlgVar))){ 
       putRcmdr("catlgVar", tclVar("catlg"))
       message("The specified design catalogue does not exist and has been replaced by the default!")
  }
  else{
    if (!"catlg" %in% class(get(tclvalue(catlgVar)))){ 
         putRcmdr("catlgVar", tclVar("catlg"))
         message("The specified design catalogue is invalid and has been replaced by the default!")
         }
  }
  designrbVariable <- tclVar(.stored.design2FrF$designrbVariable) 
  designrbFrame <- ttkframe(designFrame)
  defaultrb <- tkradiobutton(designrbFrame,text=gettext("None"),variable=designrbVariable,value="default", command=onRefresh)
  genrb <- tkradiobutton(designrbFrame,text=gettext("Specify Generators (generators option)"),variable=designrbVariable,value="gen", command=onRefresh)
  catlgrb <- tkradiobutton(designrbFrame,text=gettext("Specify catalogue name (select.catlg option)"),variable=designrbVariable,value="catlg", command=onRefresh)
  designrb <- tkradiobutton(designrbFrame,text=gettext("Specify Design name (design option)"),variable=designrbVariable,value="design", command=onRefresh)
  tkgrid(defaultrb,sticky="w")
  tkgrid(genrb,sticky="w")
  tkgrid(catlgrb,sticky="w")
  tkgrid(designrb,sticky="w")
  tkgrid(tklabel(designrbFrame,text="  "))

  genVar <- tclVar(.stored.design2FrF$genVar)
  designVar <- tclVar(.stored.design2FrF$designVar)
  genEntry <- tkentry(designFrame, textvariable=genVar, width=50)
  
  tkgrid(designrbFrame,columnspan=4,sticky="n")
  
  ## prepare list of catalogue names
  ## catlgliste (without s) is list of catlg entries of current catalogue
  putRcmdr("catlgsliste", listCatlgs())
  catlgEntry <- ttkcombobox(designFrame, textvariable=catlgVar, width=24, values=catlgsliste)
    tkbind(catlgEntry, "<<ComboboxSelected>>", catlgnamechange)
  
  
  ## prepare list of design names, based on current setting for the catalogue
        pick <- nfac(get(tclvalue(catlgVar)))==as.numeric(tclvalue(nfacVar))
        if (tclvalue(resVar)=="IV") pick <- pick & res(get(tclvalue(catlgVar))) > 3 else
        if (tclvalue(resVar)=="V+") pick <- pick & res(get(tclvalue(catlgVar))) > 4
        if (tclvalue(nrunEntryVariable)=="1") pick <- pick & nruns(get(tclvalue(catlgVar))) == as.numeric(tclvalue(nrunVar))
        if (any(pick)) putRcmdr("catlgliste", names(get(tclvalue(catlgVar)))[pick])
        else putRcmdr("catlgliste", "")
        if (!tclvalue(designVar) %in% catlgliste & length(catlgliste) > 0) designVar <- tclVar(catlgliste[1])
  
  
  designEntry <- ttkcombobox(designFrame, textvariable=designVar, values=catlgliste, width=24, state="readonly")
  
  if (tclvalue(designrbVariable)=="gen") 
  {   tkgrid(tklabel(designFrame,text=gettext("Type in generators")), columnspan=2, sticky="w")
      tkgrid(genEntry, columnspan=2)
      tkgrid(tklabel(designFrame,text=gettext("comma-separated column numbers of Yates matrix (e.g. 7, 13, 11)\nor comma-separated interaction columns (e.g. ABC, ACD, ABD)"), wraplength=500),columnspan=2)
  }
  if (tclvalue(designrbVariable)=="catlg") 
  {  tkgrid(catlab,catlgEntry,sticky="nw")
  }
  if (tclvalue(designrbVariable)=="design") 
  {  tkgrid(catlab,catlgEntry,sticky="nw")
     tkgrid(tklabel(designFrame,text=gettext("Select design\n(preselected according to numbers of runs and factors and resolution):"), wraplength=500, justify="left"),columnspan=2, sticky="w")
      tkgrid(designEntry, sticky="w", padx=15)
  }

## preparations for bottom frame
bottomFrame <- tkframe(topdes2)
specialcbVariable <- tclVar(.stored.design2FrF$cbInitials[7])
specialcb <- ttkcheckbutton(bottomFrame,text=gettext("Activate Special Choices"),
    variable=specialcbVariable, command = onSpecialcb)


helptab1Button <- buttonRcmdr(tab1, text = gettext("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab1, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab1Button, takefocus=0)


### Finalize tab1
tkgrid(tklabel(tab1, text="Name of new design"), nameEntry, helptab1Button, sticky="w", pady=10)
tkgrid.configure(helptab1Button, sticky="e")

## grid base frame
tkgrid(nrunlab <- tklabel(baseFrame, text=gettext("Number of runs")), nrunEntry, nruncb, sticky="w")
## omitted nfaccb, on form, nfactors must always be specified
tkgrid(nfaclab <- tklabel(baseFrame, text=gettext("Number of factors")), nfacEntry, sticky="w")
tkgrid(ncenterlab <- tklabel(baseFrame, text=gettext("Number of center points")), ncenterEntry, sticky="w")
tkgrid.configure(ncenterlab, pady=10)
tkgrid(tklabel(baseFrame, text=gettext("Number of blocks")), nblockEntry, aliasblock2fi, sticky="w")
tkgrid(tklabel(baseFrame, text=gettext("Replications")), nrepEntry, repeat.onlycb, sticky="w",pady=5)
tkgrid(tklabel(baseFrame, text="You normally do not need to change randomization settings"),sticky="w",columnspan=3)
tkgrid(seedlab <- tklabel(baseFrame, text=gettext("Seed for randomization")), seedEntry, 
       randomizecb, sticky="w")



## base FrF2
## expert choices for FrF2 in case of special
    if (as.logical(as.numeric(tclvalue(specialcbVariable)))){ 
       tkgrid(baseFrame, designFrame, sticky="w", columnspan=2)
       #tkgrid.configure(baseFrame, columnspan=2)
       } 
    else tkgrid(baseFrame,sticky="w",columnspan=2)
## design properties with lookup
tkgrid(despropFrame, sticky="w",columnspan=3, pady=10)

## Factor Details Tab
## factor details frame
### faclevelCommonVariable

## default levels frame
deflevFrame <- ttklabelframe(tab2,text="Default levels")
faclevelCommonVariable <- tclVar(.stored.design2FrF$cbInitials[4])
faclevelCommonButton <- ttkcheckbutton(deflevFrame,text=gettext("Common factor levels"),
    variable=faclevelCommonVariable,command=onRefresh)
tkconfigure(faclevelCommonButton,takefocus=0)
putRcmdr("level1Var", tclVar(.stored.design2FrF$level1Var))
level1Entry <- ttkentry(deflevFrame, width="20", textvariable=level1Var)
    tkconfigure(level1Entry,takefocus=0)
    tkbind(level1Entry, "<FocusIn>", level1enter)
    tkbind(level1Entry, "<FocusOut>", level1change)

putRcmdr("level2Var", tclVar(.stored.design2FrF$level2Var))
level2Entry <- tkentry(deflevFrame, width="20", textvariable=level2Var)
    tkconfigure(level2Entry,takefocus=0)
    tkbind(level2Entry, "<FocusIn>", level2enter)
    tkbind(level2Entry, "<FocusOut>", level2change)
tkgrid(faclevelCommonButton,sticky="w",columnspan=3, padx=15)
faclevCommonLab<-tklabel(deflevFrame,text=gettext("CAUTION: Checking this box overwrites all custom factor levels."))
if (!as.logical(as.numeric(tclvalue(faclevelCommonVariable)))) tkgrid(faclevCommonLab,sticky="w", columnspan=3,pady=5)
tkgrid(tklabel(deflevFrame, text=gettext("First Level")),tklabel(deflevFrame, text=gettext("Second Level")),sticky="e",padx=15)
tkgrid(level1Entry, level2Entry, sticky="e",padx=15)

## factor details
## values as vectors
facnamlistt <- .stored.design2FrF$facnamlist
if (as.logical(as.numeric(tclvalue(faclevelCommonVariable)))) {
    faclev1listt <- rep(tclvalue(level1Var),tclvalue(nfacVar)) 
    faclev2listt <- rep(tclvalue(level2Var),tclvalue(nfacVar)) 
    } else{
    faclev1listt <- .stored.design2FrF$faclev1list
    faclev2listt <- .stored.design2FrF$faclev2list
    }
faclablistt <- .stored.design2FrF$faclablist
varlistshortt <- if (as.numeric(tclvalue(nfacVar))<=50) 
                 Letters[1:tclvalue(nfacVar)] else paste("F",1:tclvalue(nfacVar),sep="")

### this deletes all user entries, if user changes number of factor levels after 
### inputting everything
### and all updating facilities go wrong
#if (!length(facnamlistt)==length(varlistshortt)){
#       faclev1listt <- rep(tclvalue(level1Var),tclvalue(nfacVar))
#       faclev2listt <- rep(tclvalue(level2Var),tclvalue(nfacVar))
#       facnamlistt <- varlistshortt
#       faclablistt <- rep("",as.numeric(tclvalue(nfacVar)))
#    }
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
    enterFrame <- ttklabelframe(enterlistFrame, text=gettext("Modify factor details for selected factor"))
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
    tkgrid(tklabel(enterFrame,text=gettext("Select"),width=6),
           tklabel(enterFrame,text=gettext("Factor name"), width=15),
           tklabel(enterFrame,text=gettext("First level"), width=15),
           tklabel(enterFrame,text=gettext("Second level"), width=15),
           tklabel(enterFrame,text=gettext("Comment or label \n(for html export only)"), width=20),
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
#    tkgrid.configure(enterFrame, columnspan=5)
    ## Hoch-/Runterschieben von Einträgen ermöglichen

    downupFrame <- ttkframe(listFrame)
    moveDownButton <- buttonRcmdr(downupFrame, text = gettext("Move Down"), 
        foreground = "darkgreen", command = onDown, 
        default = "normal", borderwidth = 3, width=12)
    moveUpButton <- buttonRcmdr(downupFrame, text = gettext("Move Up"), 
        foreground = "darkgreen", command = onUp, 
        default = "normal", borderwidth = 3, width=12)
    tkgrid(moveDownButton, sticky="w")
    tkgrid(moveUpButton, sticky="w")

    tkgrid(scrollbar, facshortListBox, facnameListBox, faclev1ListBox, faclev2ListBox, faclabListBox, downupFrame, sticky = "nw")
    tkgrid.configure(scrollbar, sticky = "wns")
    tkgrid.configure(facnameListBox, sticky = "new")

    ## hard frame, to be displayed if special choices are activated
   putRcmdr("hardVar", tclVar(.stored.design2FrF$hardVar))
   hardFrame <- ttklabelframe(tab2, text=gettext("How many factors are hard to change ?"))
   hardEntry <- tkentry(hardFrame, width="3", textvariable=hardVar)
   tkconfigure(hardEntry,takefocus=0)
   
   hardlabelFrame <- ttklabelframe(hardFrame, text=gettext("WARNING"))

   tkgrid(tklabel(hardFrame,text="The first "),hardEntry,tklabel(hardFrame,text=" factors are hard to change."))
   tkgrid(tklabel(hardFrame,text="If necessary, modify the factor order on the Factor Details tab!"),columnspan=3)
   tkgrid(tklabel(hardlabelFrame,text=gettext("Only specify this, if some factors are   r e a l l y   hard to change.\nThe hard-to-change factors arrangement of experimental runs bears the risk of misjudging effects of these factors!"),wraplength=300))
   tkgrid(hardlabelFrame, sticky="w",columnspan=3)

helptab2Button <- buttonRcmdr(tab2, text = gettext("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab2, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab2Button, takefocus=0)

## finalize tab2 Factor details
    tkgrid(helptab2Button, sticky="e")
    if (as.logical(as.numeric(tclvalue(specialcbVariable)))){ 
    tkgrid(deflevFrame, hardFrame, sticky="nw")
    }
    else tkgrid(deflevFrame, sticky="nw")
    tkgrid(listFrame, columnspan=6,sticky="w",pady=10)
    tkgrid(enterlistFrame, columnspan=6,sticky="w")

## tab3
helptab3Button <- buttonRcmdr(tab3, text = gettext("Tab Help"), 
        foreground = "darkgreen", command = onHelpTabEstimable, 
        default = "normal", borderwidth = 3)
tkconfigure(helptab3Button, takefocus=0)
tkgrid(helptab3Button, sticky="e", columnspan=6)

estradioFrame <- ttklabelframe(tab3, text=gettext("Mode of requesting estimable 2-factor interactions"))
putRcmdr("estrbVariable", tclVar(.stored.design2FrF$estrbVariable))
noestrb <- tkradiobutton(estradioFrame,text=gettext("None"),variable=estrbVariable,value="none",command=onestrb)
clearrb <- tkradiobutton(estradioFrame,text=gettext("Selected interactions must be clear of aliasing with ANY 2-factor interactions"),variable=estrbVariable,value="clear",wraplength="500",justify="left",command=onestrb)
distinctrb <- tkradiobutton(estradioFrame,text=gettext("Selected interactions must be clear of aliasing with EACH OTHER"),variable=estrbVariable,value="distinct",wraplength="500",justify="left",command=onestrb)
## deactivate all fields on this tab, if this box is not checked
estlabel <- tklabel(estradioFrame, text=gettext("NOTE: Resolution entry on tab Base Settings is ignored, if choice is not None"))
tkgrid(noestrb, sticky="w")
tkgrid(clearrb, sticky="w")
tkgrid(distinctrb, sticky="w")
tkgrid(estlabel, sticky="w")

resoFrame <- ttklabelframe(tab3,text=gettext("Minimum Resolution"))
resolabel <- tklabel(resoFrame, text=gettext("Per default, at least resolution IV is required. \nThis is natural for most applications."),
                    justify="left")
res3cbVariable <- tclVar(.stored.design2FrF$cbInitials[10])
res3cb <- ttkcheckbutton(resoFrame,text=gettext("Permit a resolution III design"),variable=res3cbVariable)
tkgrid(resolabel, sticky="w", columnspan=5)
tkgrid(res3cb, sticky="w", columnspan=5)

if (tclvalue(estrbVariable)=="none") tkgrid(estradioFrame, sticky="w", columnspan=4)
else tkgrid(estradioFrame, resoFrame, sticky="w", columnspan=5)

## create empty row
tkgrid(tklabel(tab3,text="   "))
selFrame <- ttklabelframe(tab3, text=gettext("Select 2-factor interactions"))
comprradioFrame <- ttklabelframe(selFrame, text=gettext("Type of specification"))

## design from older version of RcmdrPlugin.DoE
if (is.null(.stored.design2FrF$comprrbVariable)) 
         .stored.design2FrF$comprrbVariable <- "manual"
putRcmdr("comprrbVariable", tclVar(.stored.design2FrF$comprrbVariable))
manualestrb <- tkradiobutton(comprradioFrame,text=gettext("Select manually"),
      variable=comprrbVariable,value="manual",command=oncomprestrb)
comprestrb <- tkradiobutton(comprradioFrame,text=gettext("Pre-specified structure from two groups of factors"),
      variable=comprrbVariable,value="compr",wraplength="500",justify="left",command=oncomprestrb)

varlistshortt=(strsplit(tclvalue(varlistshort)," ")[[1]])
if (tclvalue(comprrbVariable)=="manual"){
  hilf <- combn(length(varlistshortt),2)
  intaclistt <- paste(varlistshortt[hilf[1,]],varlistshortt[hilf[2,]],sep="")
  putRcmdr("est2fislist", setdiff(.stored.design2FrF$est2fislist,setdiff(.stored.design2FrF$est2fislist,intaclistt)))
    ## omit clicked already selected interactions that must be lost because of a reduction in nfactors
}
else{
  intaclistt <- varlistshortt
  putRcmdr("est2fislist", setdiff(.stored.design2FrF$est2fislist,setdiff(.stored.design2FrF$est2fislist,intaclistt)))
    ## omit clicked already selected factors that must be lost because of a reduction in nfactors
}

## design from older version of RcmdrPlugin.DoE
if (is.null(.stored.design2FrF$comprclassVar)) 
         .stored.design2FrF$comprclassVar <- "3: all interactions of group 1"
putRcmdr("comprclassVar", tclVar(.stored.design2FrF$comprclassVar))

#resEntry <- tkentry(despropframe, textvariable=resVar)
    putRcmdr("comprclassEntry", ttkcombobox(comprradioFrame, textvariable=comprclassVar, width=50, 
         values=c("1: interactions within group1",
                  "2: interactions within groups 1 and groups 2",
                  "3: all interactions of group 1", 
                  "4: interactions between groups 1 and 2"), state="readonly"))
    #tkbind(comprclassEntry, "<<ComboboxSelected>>", onRefresh)

tkgrid(manualestrb, sticky="w")
tkgrid(comprestrb, comprclassEntry, sticky="w")
tkgrid(comprradioFrame, sticky="w", columnspan=6)
if (!tclvalue(comprrbVariable)=="compr") tkconfigure(comprclassEntry, state="disabled")

estbuttonFrame <- ttkframe(selFrame)
selectButton <- buttonRcmdr(estbuttonFrame, text = gettext(">"), 
        foreground = "darkgreen", command = onSelect, 
        default = "normal", borderwidth = 3)
tkgrid(selectButton)
deselectButton <- buttonRcmdr(estbuttonFrame, text = gettext("<"), 
        foreground = "darkgreen", command = onDeselect, 
        default = "normal", borderwidth = 3)
tkgrid(deselectButton)
onestrb.worefresh()

putRcmdr("notest2fislist", setdiff(intaclistt, est2fislist))

## define variable list boxes for selection
## intaclistt is the master list

## make sure that both lists are of equal length and long enough
## ## make selection offering depend on choice of compromise settings
     if (tclvalue(comprrbVariable)=="compr"){
        TITEL.LINKS <- "Group 1 (at least 1 element)"
        TITEL.RECHTS <- "Group 2 (at least 1 element)"
     }
     else {
        TITEL.LINKS <- "Available 2-factor interactions"
        TITEL.RECHTS <- "Selected 2-factor interactions"
     }
putRcmdr("est2fis", variableListBox(selFrame, variableList=intaclistt, listHeight=15, 
    title=TITEL.RECHTS, selectmode="multiple"))
putRcmdr("notest2fis", variableListBox(selFrame, variableList=intaclistt, listHeight=15, 
    title=TITEL.LINKS, selectmode="multiple"))
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist


maxtimeFrame <- ttklabelframe(selFrame, text=gettext("Limit search time"))
maxtimelabel <- tklabel(maxtimeFrame, 
    text=gettext("The search can take very long.\nIf it times out unsuccessfully, \nexpert users may try the command line mode of function FrF2 that allows more control options."),
    justify="left", wraplength="280")
maxtimeentrylabel <- tklabel(maxtimeFrame, text=gettext("Maximum search time in seconds"))
maxtimeVar <- tclVar(.stored.design2FrF$maxtimeVar)
maxtimeEntry <- tkentry(maxtimeFrame, textvariable=maxtimeVar)
tkgrid(maxtimeentrylabel, sticky="w")
tkgrid(maxtimeEntry, sticky="e")
tkgrid(maxtimelabel,sticky="w")
## create empty row
tkgrid(tklabel(tab3,text="   "))
#tkgrid(maxtimeFrame, sticky="w", columnspan=5)

if (tclvalue(estrbVariable)=="distinct") {tkgrid(notest2fis$frame, estbuttonFrame, est2fis$frame, maxtimeFrame, sticky="w")
     tkgrid.configure(maxtimeFrame, sticky="e")
}
else tkgrid(notest2fis$frame, estbuttonFrame, est2fis$frame, sticky="w")
tkgrid(selFrame, sticky="ew", columnspan=6)

#if (tclvalue(estrbVariable)=="none"){
#    tkconfigure(selectButton, state="disabled")
#    tkconfigure(deselectButton, state="disabled")
#    tkconfigure(maxtimeentrylabel, state="disabled")
#    tkconfigure(maxtimeEntry, state="disabled")
#    tkconfigure(maxtimelabel, state="disabled")
#    tkconfigure(resolabel, state="disabled")
#    tkconfigure(res3cb, state="disabled")
#}
#else{
#    tkconfigure(selectButton, state="normal")
#    tkconfigure(deselectButton, state="normal")
#    tkconfigure(resolabel, state="normal")
#    tkconfigure(res3cb, state="normal")
#    if (tclvalue(estrbVariable)=="distinct") {
#    tkconfigure(maxtimeentrylabel, state="normal")
#    tkconfigure(maxtimeEntry, state="normal")
#    tkconfigure(maxtimelabel, state="normal")
#        } 
#    else {
#    tkconfigure(maxtimeentrylabel, state="disabled")
#    tkconfigure(maxtimeEntry, state="disabled")
#    tkconfigure(maxtimelabel, state="disabled")
#        }
#}

tkgrid(tklabel(tab4, text=gettext("This tab will accomodate block or split-plot functionality.")), sticky="n")

### tab5
## define frames for special activities on specials tab
specialFrame <- ttklabelframe(tab5,text=gettext("Any special non-standard requirements ?"))

### widgets for special frame
putRcmdr("specialrbVariable", tclVar(.stored.design2FrF$specialrbVariable) )
## for some reason, doesn't work with initialValue directly as text string
  specialrbFrame <- ttkframe(specialFrame)
  nonerb <- tkradiobutton(specialrbFrame,text=gettext("None"),variable=specialrbVariable,value="none", command=onRefresh)
  hardrb <- tkradiobutton(specialrbFrame,text=gettext("(Some) Hard to change factor(s)"),variable=specialrbVariable,value="hard", command=onRefresh)
  debarrb <- tkradiobutton(specialrbFrame,text=gettext("Some combinations impossible"),variable=specialrbVariable,value="debar", command=onRefresh)
  tkgrid(nonerb,sticky="w")
  tkgrid(hardrb,sticky="w")
  tkgrid(debarrb,sticky="w")
  tkgrid(tklabel(specialrbFrame,text="  "))

tkgrid(specialrbFrame, columnspan=4,sticky="w")
if (tclvalue(specialrbVariable)=="none") tkgrid(specialFrame, sticky="nw",columnspan=5)

#if (tclvalue(specialrbVariable)=="hard"){
#   tkgrid(specialFrame, hardlabelFrame, sticky="w",columnspan=5)
#   tkgrid(hardFrame, sticky="w")
#}

if (tclvalue(specialrbVariable)=="debar"){
   debarlabelFrame <- ttklabelframe(tab5,text=gettext("WARNING"))
   
   nrunEntryDebar <- ttkentry(tab5, textvariable=nrunVar)
   debarFrame <- ttklabelframe(tab5, text=gettext("Conditions for excluding combinations"))
   tkgrid(tklabel(debarFrame,text="must think of a simple way of specifying conditions with clicking them together (or typing at users choice)",wraplength=500))
   tkgrid(tklabel(debarlabelFrame,text="If some experimental combinations are impossible, standard designs cannot be used. Please try to avoid this.\n\nIf it is unavoidable, estimable effects must be specified on the Estimable Model tab, with the third radio button in the top left frame activated.", wraplength=500, justify="left"),sticky="w") 
   tkgrid(specialFrame, debarlabelFrame, sticky="w")

   tkgrid(ttklabel(tab5,text="  "),sticky="w")
   tkgrid(ttklabel(tab5,text=gettext("Number of runs: ")),nrunEntryDebar,sticky="w")
   tkgrid(ttklabel(tab5,text=gettext("The number of runs need not be a power of 2 any more, as the design is determined as a non-orthogonal D-optimal design."), justify="left", wraplength=500),sticky="w",columnspan=2)
   tkgrid(ttklabel(tab5,text="  "),sticky="w")
      
   tkgrid(debarFrame, columnspan=2,sticky="w")
   ### here, we need the possibility to enter conditions to debar combinations
}

## finalize tab5


## tab6 for exporting
helptab6Button <- buttonRcmdr(tab6, text = gettext("Tab Help"), 
        foreground = "darkgreen", command = onHelpTab6, 
        default = "normal", borderwidth = 3)

exportlabVar <- nameVar
exportlab <- ttklabel(tab6, textvariable=exportlabVar)
tkgrid(ttklabel(tab6,text="Current design to be saved:"),exportlab,helptab6Button,sticky="w", pady=10) 

## radio buttons for choosing export type
etradioFrame <- ttklabelframe(tab6, text=gettext("(How to) Export ?"))
etyperbVariable <- tclVar(.stored.design2FrF$etyperbVariable)
noexprb <- tkradiobutton(etradioFrame,text=gettext("no export"),variable=etyperbVariable,value="none")
allrb <- tkradiobutton(etradioFrame,text=gettext("all file types"),variable=etyperbVariable,value="all")
rdarb <- tkradiobutton(etradioFrame,text=gettext("rda only"),variable=etyperbVariable,value="rda")
htmlrb <- tkradiobutton(etradioFrame,text=gettext("html and rda"),variable=etyperbVariable,value="html")
csvrb <- tkradiobutton(etradioFrame,text=gettext("csv and rda"),variable=etyperbVariable,value="csv")
tkgrid(noexprb, sticky="w")
tkgrid(allrb, sticky="w")
tkgrid(rdarb, sticky="w")
tkgrid(htmlrb, sticky="w")
tkgrid(csvrb, sticky="w")

## radio buttons for choosing export decimal separator
decimalradioFrame <- ttklabelframe(tab6, text=gettext("Decimal Separator ?"))
decimalrbVariable <- tclVar(.stored.design2FrF$decimalrbVariable)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettext("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettext("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettext(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(tab6, text=gettext("Storage Directory"))
putRcmdr("dirVar", tclVar(.stored.design2FrF$dirVar))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettext("Change directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, dirButton, sticky="w")
tkgrid.configure(dirEntry,padx=15)

## export file name
putRcmdr("fileVar", tclVar(.stored.design2FrF$fileVar))
fileEntry <- tkentry(tab6, width="20", textvariable=fileVar)
efnamelabel <- tklabel(tab6,text=gettext("Export file names: name below with appropriate endings (html or csv, and rda)"))
putRcmdr("replacecbVariable", tclVar(.stored.design2FrF$cbInitials[8]))
replacecb <- ttkcheckbutton(tab6,text=gettext("Replace file(s), if exists"),variable=replacecbVariable)

## always grid details, as otherwise default file name does not work
## design name info and help button have already been gridded above
tkgrid(etradioFrame, decimalradioFrame, sticky="w")
tkgrid(dirFrame, sticky="w", columnspan=5)
tkgrid.configure(dirFrame, pady=15)
tkgrid(efnamelabel, sticky="w", columnspan=5)
tkgrid(fileEntry, sticky="w", columnspan=5)
tkgrid(replacecb, sticky="w", columnspan=5)


## add buttons outside the notebook
buttonFrame <- tkframe(topdes2)
## die sind aber nicht dunkelgruen ...
refreshButton <- buttonRcmdr(buttonFrame, text = gettext("Refresh form"), 
        foreground = "darkgreen", width = "12", command = onRefresh, 
        default = "normal", borderwidth = 3)
storeButton <- buttonRcmdr(buttonFrame, text = gettext("Store form"), 
        foreground = "darkgreen", width = "12", command = onStore, 
        default = "normal", borderwidth = 3)
loadButton <- buttonRcmdr(buttonFrame, text = gettext("Load form"), 
        foreground = "darkgreen", width = "12", command = onLoad, 
        default = "normal", borderwidth = 3)
resetButton <- buttonRcmdr(buttonFrame, text = gettext("Reset form"), 
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

OKCancelHelp(window=topdes2, helpSubject="Menu.FrF2level")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

#tkbind(specialcb, "<ButtonRelease-1>", onReset)
tkgrid(specialcb, sticky="e")
tkconfigure(specialcb, takefocus=0)

tkgrid(buttonsFrame, bottomFrame, sticky="ew")


### relations among widgets
if (!as.logical(as.numeric(tclvalue(randomizeVariable)))){
tkconfigure(seedEntry, state="disabled")
tkconfigure(seedlab, state="disabled")
}else {
tkconfigure(seedEntry, state="normal")
tkconfigure(seedlab, state="normal")
}
if (!as.logical(as.numeric(tclvalue(nrunEntryVariable)))){
tkconfigure(nrunEntry, state="disabled")
tkconfigure(nrunlab, state="disabled")
}else {
tkconfigure(nrunEntry, state="normal")
tkconfigure(nrunlab, state="normal")
}
##if (!as.logical(as.numeric(tclvalue(nfacEntryVariable)))){
##tkconfigure(nfacEntry, state="disabled")
##tkconfigure(nfaclab, state="disabled")
##}else {
##tkconfigure(nfacEntry, state="normal")
##tkconfigure(nfaclab, state="normal")
##}

           
if (as.logical(as.numeric(tclvalue(specialcbVariable)))){
  tcl(tn, "tab", 2, state="normal")
  tcl(tn, "tab", 3, state="hidden")
  tcl(tn, "tab", 4, state="hidden")
  ## hide tabs that are not yet usable
}
else {
  tcl(tn, "tab", 2, state="hidden")
  tcl(tn, "tab", 3, state="hidden")
  tcl(tn, "tab", 4, state="hidden")
}


#### former tab5: special
## do last in order to have it overwrite the others rather than 
## have others disable special again
#if (as.logical(as.numeric(tclvalue(specialcbVariable)))){
#if (tclvalue(specialrbVariable) == "block" || tclvalue(specialrbVariable) == "splitplot")
#RcmdrTkmessageBox(gettext("You cannot declare impossible combinations for blocked or split-plot designs. Please resolve."), 
#     icon="error", type="ok", title="Conflict")
#  tcl(tn, "tab", 4, state="normal")
#if (!tclvalue(hardVar) == "0")
#     RcmdrTkmessageBox(gettext("You cannot declare impossible combinations for designs with hard-to-change factors. Please resolve."), 
#     icon="error", type="ok", title="Conflict")
#  tcl(tn, "tab", 4, state="normal")
#}else {
#  tcl(tn, "tab", 4, state="hidden")
#}

#if (as.numeric(tclvalue(hardVar))>0){
#if (tclvalue(specialrbVariable) == "block" || tclvalue(specialrbVariable) == "splitplot")
#RcmdrTkmessageBox(gettext("You cannot combine hard-to-change variables with blocked or split-plot designs. (A split-plot design in itself can often take care of hard-to-change variables.) Please resolve."), 
#     icon="error", type="ok", title="Conflict")
#if (tclvalue(specialrbVariable) == "estimable")
#RcmdrTkmessageBox(gettext("You cannot combine hard-to-change variables with estimable models. Please resolve."), 
#     icon="error", type="ok", title="Conflict")
#}
if (exists("activestab.tn", where="RcmdrEnv")){
                tcl(tn, "select", activestab.tn)
                rm(activestab.tn, pos="RcmdrEnv")
                }

dialogSuffix(window=topdes2, rows=2, columns=2, focus=tn, bindReturn=FALSE)

}
###
# End of Menu.FrF2level
###
