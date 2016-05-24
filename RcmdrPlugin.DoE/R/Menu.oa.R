## two instances of assign replaced by justDoIt

Menu.oa <- function(){
initializeDialogDoE(title=gettextRcmdr("Create main effects design from orthogonal array ..."))   
     ## function initializeDialogDoE assumes topdes2 as windowname
     ## last stored top left corner for window is stored under topleft2xy
     ## onRefresh still makes window walk a little
     

if (exists("curindex", where="RcmdrEnv")) rm(curindex, pos="RcmdrEnv")

if (!exists(".stored.designoa",where="RcmdrEnv")) 
           assign(".stored.designoa", .default.designoa,pos="RcmdrEnv")
           ## nameVar, nrunVar, nfacVar, nrepVar
           ## cbInitials containing repeat.onlyVariable, randomizeVariable, 
           ##                       facnamesAutoVariable, faclevelsCommonVariable, 
           ##                       nrunEntryVariable, usenlevelscbVariable
           ##                       specialcbVariable, replacecbVariable, parentcbVariable,
           ##                       res3cbVariable
           ## level1Var, level2Var, seedVar, specialrbVariable, hardVar, genVar, 
           ## catlgVar, designVar, designrbVariable, destyperbVariable
           ## resVar, qualcritrbVariable, facnamlist,nlevlist,faclevlist, faclablist
           ## etyperbVariable, decimalrbVariable, dirVar, fileVar
           ## fromVar, toVar

putRcmdr("nrunslist",c("NULL",as.character(unique(oacat$nruns))))
idlist <- as.character(oacat$name)
idlist <- idlist[-which(oacat$lineage=="" & 
              oacat$nruns==sapply(strsplit(idlist,".", fixed=TRUE), function(obj) 
              prod(apply(matrix(as.numeric(obj[-1]), nrow=2), 2, function(obj2) obj2[1]^obj2[2]))))]
idlist <- c("NULL",idlist)
putRcmdr("idlist", idlist)

## define called functions
 infoClose <- function(){
     putRcmdr("infotxt",tclVar(""))
 }
 
 onHelpTab1 <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Menu.oaTab1"))
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

onInspect <- function(){
     if (tclvalue(fromVar)=="" & tclvalue(toVar)=="") {
        if (!tclvalue(nrunVar)=="NULL") nrunshow <- paste("c(",tclvalue(nrunVar),", 100000),")
        else nrunshow <- "\"all\", "
     }
     else{
        if (tclvalue(fromVar) == "") nrunshow <- paste("c(0,", tclvalue(toVar),"),")
        else {if (tclvalue(toVar) =="") nrunshow <- paste("c(",tclvalue(fromVar),", 100000),")
              else if (tclvalue(toVar)==tclvalue(fromVar)) nrunshow <- paste(tclvalue(toVar),",", sep="")
                   else nrunshow <- paste("c(",tclvalue(fromVar),",",tclvalue(toVar),"),")
             }
     }
     if (as.logical(as.numeric(tclvalue(usenlevelscbVariable))))
         nlevelsshow <- paste("nlevels = c(",paste(as.character(tclObj(nlevlist)),collapse=","),"),")
         else nlevelsshow <- ""
     command <- paste("show.oas(nruns = ", nrunshow, nlevelsshow, 
        " parents.only = ", as.logical(as.numeric(tclvalue(parentcbVariable))), 
        ", show = Inf)")
     hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            logger(gettextRcmdr("Something went wrong, perhaps invalid values for number of runs?"))
            return()
            }
      ## logger(command)  ## braucht man nur für justDoIt
      tkgrab.release(topdes2)
      doItAndPrint(command)
      }

storeRcmdr <- function(){
        hilf <- list(nameVar=tclvalue(nameVar),idVar = tclvalue(idVar), 
        nrunVar=tclvalue(nrunVar),nfacVar=tclvalue(nfacVar),nrepVar=tclvalue(nrepVar), 
        minrdfVar=tclvalue(minrdfVar),
        optimVar = tclvalue(optimVar),
        cbInitials = c(tclvalue(repeat.onlyVariable), tclvalue(randomizeVariable),
                       tclvalue(colnospecifyVariable),0,
                       1,tclvalue(usenlevelscbVariable),
                       0,tclvalue(replacecbVariable),tclvalue(parentcbVariable),
                       0
                       ),
        seedVar=tclvalue(seedVar),
        facnamlist=as.character(tclObj(facnamlist)),
        nlevlist=as.character(tclObj(nlevlist)),
        colnolist=as.character(tclObj(colnolist)),
        faclevlist=as.character(tclObj(faclevlist)),
        faclablist=as.character(tclObj(faclablist)),
        etyperbVariable=tclvalue(etyperbVariable),
        decimalrbVariable=tclvalue(decimalrbVariable),
        dirVar=tclvalue(dirVar), fileVar=tclvalue(fileVar),
        fromVar=tclvalue(fromVar), toVar=tclvalue(toVar))
        class(hilf) <- c("menu.designoa","list")
        putRcmdr(".stored.designoa",hilf)
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
            errorCondition(window=topdes2,recall=Menu.oa, 
                    message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Object"))))
            {
              errorCondition(window=topdes2,recall=Menu.oa, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
        if (any(getRcmdr(".stored.designoa")$nlevlist=="") )
          {
              errorCondition(window=topdes2,recall=Menu.oa, 
              gettextRcmdr("Factor details are not completely specified."))
              return()
          }
        if (any(getRcmdr(".stored.designoa")$faclevlist=="") )
          {
              errorCondition(window=topdes2,recall=Menu.oa, 
              gettextRcmdr("Factor details are not completely specified."))
              return()
          }
    ###  further error messages with return to menu ?

   
    textfactornameslist.forcommand <- paste("factor.names=list(",paste(paste(as.character(tclObj(facnamlist)),"=c(",
                            sapply(strsplit(as.character(tclObj(faclevlist))," "),function(obj) paste(dquote(obj),collapse=",")), 
                            ")",sep=""),collapse=","),")")

    
    columns <- dquote("order")
    if (as.logical(as.numeric(as.character(tclvalue(colnospecifyVariable)))))
       columns <- paste("c(",paste(as.character(tclObj(colnolist)), collapse=","),")")
    
    if (columns == dquote("order") & !tclvalue(optimVar)=="none")
       columns <- dquote(tclvalue(optimVar))
    
    command <- paste("oa.design(ID=",tclvalue(idVar),",nruns=",tclvalue(nrunVar),
                     ",nfactors=",tclvalue(nfacVar),",replications=",
                  tclvalue(nrepVar),",repeat.only=",as.logical(as.numeric(tclvalue(repeat.onlyVariable))),
                  ",randomize=",as.logical(as.numeric(tclvalue(randomizeVariable))),",seed=",tclvalue(seedVar),
                  ",nlevels=c(", paste(as.character(tclObj(nlevlist)),collapse=","),
                  "),",textfactornameslist.forcommand, ", columns =", columns,
                  ", min.residual.df=", tclvalue(minrdfVar), ")") 

        logger(gettextRcmdr("## Trying to create design ... "))
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=topdes2,recall=Menu.oa, message=gettextRcmdr(hilf))
             return()
            }
        logger(paste(name, "<-", command))
        logger("## creator element of design.info will be different, when using the command line command!")
        ## change creator to contain menu settings
        hilfatt <- design.info(hilf)
        hilfatt$creator <- .stored.designoa
        class(hilfatt$creator) <- c("menu.designoa", "list")
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
            errorCondition(window=topdes2,recall=Menu.oa, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        }
        rm(activestab.tn, pos="RcmdrEnv")
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
  }

listDesignoa <- function (envir = .GlobalEnv, ...) 
{
    Vars <- ls(envir = envir, all.names = TRUE)
    Vars[which(sapply(Vars, function(.x){
               aus <- FALSE
               if ("menu.designoa" %in% class(get(.x, envir = envir))) aus <- TRUE
               else if ("design" %in% class(get(.x, envir = envir)))
                    if ("menu.designoa" %in% class(design.info(get(.x, envir = envir))$creator))
                       aus <- TRUE
               aus
               }))]
}


onLoad <- function(){
    ## seems to work now, needs to be tested!
        hilf <- listDesignoa()
        if (length(hilf)==0) {
            tkmessageBox(message=gettextRcmdr("There are no stored design inputs in this session."),
               icon="error", type="ok", title="no stored design inputs")
            return()
            }
    putRcmdr("deschooseoa",tktoplevel())
    tkwm.title(deschooseoa, gettextRcmdr("Choose stored design form"))
    position <- if (is.SciViews()) 
        -1
    else position <- "+50+50"
    tkwm.geometry(deschooseoa, position)
    putRcmdr("lb", variableListBox(deschooseoa, variableList=hilf, title="Choose stored design form"))
        tkgrid(lb$frame)
    onOK <- function() {
        putRcmdr(".stored.designoa",get(lb$varlist[as.numeric(tclvalue(tcl(lb$listbox, "curselection")))+1]))
        if ("design" %in% class(getRcmdr(".stored.designoa"))) 
            putRcmdr(".stored.designoa", design.info(getRcmdr(".stored.designoa"))$creator)
        tkfocus(CommanderWindow())
        tkdestroy(topdes2)
        tkdestroy(deschooseoa)
        Menu.oa()
    }
    OKCancelHelp(window=deschooseoa)
    tkgrid(buttonsFrame, sticky="s")
    dialogSuffix(window=deschooseoa, rows=1, columns=1, 
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
          Menu.oa()
}

onStore <- function(){
        ## Speichernamen abfragen und hier ermöglichen (statt stored.designoa)
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
        ## replace assign by justDoIt; assign(savename.RcmdrPlugin.DoE, getRcmdr(".stored.designoa"), envir=.GlobalEnv)
        justDoIt(paste(savename.RcmdrPlugin.DoE, "<- getRcmdr(\".stored.designoa\")"))
        message(gettextRcmdr("inputs have been stored"))
        }
}

onReset <- function(){
        assign(".stored.designoa",.default.designoa,pos="RcmdrEnv")
        tkfocus(CommanderWindow())
  tkdestroy(topdes2)
  Menu.oa()
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
           putRcmdr("colnolist", tclVar(as.character(tclObj(colnolist))[1:nfacnew]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[1:nfacnew]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[1:nfacnew]))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(nlevListBox, listvariable=nlevlist, height=min(10,nfacnew))
           tkconfigure(colnoListBox, listvariable=colnolist, height=min(10,nfacnew))
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
           putRcmdr("colnolist", tclVar(c(as.character(tclObj(colnolist)),
                               rep("",nfacnew-nfacold))))
           putRcmdr("faclevlist", tclVar(c(as.character(tclObj(faclevlist)),
                               rep("",nfacnew-nfacold))))
           putRcmdr("faclablist", tclVar(c(as.character(tclObj(faclablist)),
                               rep("",nfacnew-nfacold))))
           tkconfigure(facshortListBox, listvariable=varlistshort, height=min(10,nfacnew))
           tkconfigure(fsel, values=varlistshortt)   
           tkconfigure(facnameListBox, listvariable=facnamlist, height=min(10,nfacnew))
           tkconfigure(nlevListBox, listvariable=nlevlist, height=min(10,nfacnew))
           tkconfigure(colnoListBox, listvariable=colnolist, height=min(10,nfacnew))
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
        putRcmdr("curcolno", tclVar(as.character(tclObj(colnolist))[selpos]))
        putRcmdr("curflev", tclVar(as.character(tclObj(faclevlist))[selpos]))
        putRcmdr("curflab", tclVar(as.character(tclObj(faclablist))[selpos]))
        tkconfigure(fnam, textvariable=curfnam)
        tkconfigure(nlev, textvariable=curnlev)
        tkconfigure(colno, textvariable=curcolno)
        tkconfigure(flev, textvariable=curflev)
        tkconfigure(flab, textvariable=curflab)
    }
    
    idsel<-function(){
        #### aendert die in der Textbox dargestellte Auswahl
        #### ruiniert aber leider auch wieder die korrekte Ueberschreibung der Werte
        putRcmdr("idpos", as.numeric(tclvalue(tcl(idEntry, "current"))) + 1)
        if (idpos==1) {
          putRcmdr("colnospecifyVariable", tclVar("0"))
          putRcmdr("colnolist", tclVar(rep("",as.numeric(tclvalue(nfacVar)))))
          tkconfigure(colnoListBox, listvariable=colnolist)
          tkconfigure(colnospecifycb, variable=colnospecifyVariable)
        }
        else {
          hilf <- eval(parse(text=paste("oa.design(",tclvalue(idVar),",randomize=FALSE)")))
          putRcmdr("nrunVar", tclVar(nrow(hilf)))
          if (as.numeric(tclvalue(nfacVar)) > ncol(hilf)) 
              tkmessageBox(message=gettextRcmdr("The chosen design cannot accomodate the chosen number of factors!"), 
               icon="error", type="ok", title=gettextRcmdr("Too many factors requested for this design"))
               }
        onRefresh()
    }
    
    nrunsel<-function(){
        #### aendert die in der Textbox dargestellte Auswahl
        #### ruiniert aber leider auch wieder die korrekte Ueberschreibung der Werte
        putRcmdr("nrunpos", as.numeric(tclvalue(tcl(nrunEntry, "current")))+1)
        if (!tclvalue(idVar)=="NULL") 
           if (nrunpos>1 & !nrow(eval(parse(text=paste("oa.design(",tclvalue(idVar),",randomize=FALSE)"))))>=tclvalue(nrunVar)) 
              tkmessageBox(message=gettextRcmdr("mismatch between chosen design and number of runs!"), 
               icon="error", type="ok", title=gettextRcmdr("Number of runs does not match chosen design"))
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
        else tkmessageBox(message=gettextRcmdr("invalid name!"),
             icon="error", type="ok", title=gettextRcmdr("Invalid factor name"))
    }
    colnochange <- function(){
        ## selpos known from factorsel
        if (length(as.character(tclObj(curcolno)))==1){
          hilf <- as.character(tclObj(colnolist))
          hilf[selpos] <- tclvalue(curcolno)
          putRcmdr("colnolist",tclVar(hilf))
          tkconfigure(colnoListBox, listvariable=colnolist)
          nncol <- ncol(eval(parse(text=tclvalue(idVar))))
         if (is.na(as.numeric(tclvalue(curcolno))) | !floor(as.numeric(tclvalue(curcolno)))==as.numeric(tclvalue(curcolno))
            | (as.numeric(tclvalue(curcolno)) > nncol )) 
            tkmessageBox(message=gettextRcmdr("an integer number from 1 to ", nncol, " was expected, please correct!"),
            icon="error", type="ok", title=gettextRcmdr("Invalid column number"))
        else{
          hilf <- as.numeric(as.character(tclObj(curcolno)))
          nnlev <- length(unique(eval(parse(text=tclvalue(idVar)))[,hilf]))
          putRcmdr("curnlev", tclVar(nnlev))
          tkconfigure(nlev, textvariable=curnlev)
          nlevchange()}
        }
        else tkmessageBox(message=gettextRcmdr("Empty entries or entries with blanks are not permitted, please correct!"),
            icon="error", type="ok", title=gettextRcmdr("Invalid column number"))
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
            icon="error", type="ok", title="Invalid number of levels")
        if (is.na(as.numeric(tclvalue(curnlev))) | !floor(as.numeric(tclvalue(curnlev)))==as.numeric(tclvalue(curnlev))) 
            tkmessageBox(message=gettextRcmdr("an integer number was expected, please correct!"),
            icon="error", type="ok", title=gettextRcmdr("Invalid number of levels"))
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
          else tkmessageBox(message=paste(gettextRcmdr(nlev," levels entries needed, \nyou specified ", nlevspec, " levels, please correct!"),collapse=" "),
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
    checkIndexColno <- function(){
        putRcmdr("curindex", as.numeric(tcl(colnoListBox,"curselection"))+1)
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
           putRcmdr("colnolist", tclVar(as.character(tclObj(colnolist))[orderUp]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[orderUp]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderUp]))
           tkconfigure(nlevListBox, listvariable=nlevlist)
           tkconfigure(colnoListBox, listvariable=colnolist)
           tkconfigure(faclevListBox, listvariable=faclevlist)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex-1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(nlevListBox,"selection","set",curindex-1)
           tcl(colnoListBox,"selection","set",curindex-1)
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
           putRcmdr("colnolist", tclVar(as.character(tclObj(colnolist))[orderDown]))
           putRcmdr("faclevlist", tclVar(as.character(tclObj(faclevlist))[orderDown]))
           putRcmdr("faclablist", tclVar(as.character(tclObj(faclablist))[orderDown]))
           tkconfigure(nlevListBox, listvariable=nlevlist)
           tkconfigure(colnoListBox, listvariable=colnolist)
           tkconfigure(faclevListBox, listvariable=faclevlist)
           tkconfigure(faclabListBox, listvariable=faclablist)
           tkconfigure(facnameListBox, listvariable=facnamlist)
           putRcmdr("curindex", curindex+1)
           indexchange()
           tcl(facshortListBox,"selection","set",curindex-1)
           tcl(nlevListBox,"selection","set",curindex-1)
           tcl(colnoListBox,"selection","set",curindex-1)
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
putRcmdr("nameVar", tclVar(.stored.designoa$nameVar))
nameEntry <- tkentry(nameFrame, width="20", textvariable=nameVar)
    tkbind(nameEntry, "<FocusIn>", nameenter)
    tkbind(nameEntry, "<FocusOut>", namechange)

nfacVar <- tclVar(.stored.designoa$nfacVar)
nfacEntry <- tkentry(baseFrame, width="8", textvariable=nfacVar)
nfacHint <- ttklabel(baseFrame, text="(required!)", foreground="#888888")
tkbind(nfacEntry,"<FocusOut>",nfacchange)

idVar <- tclVar(.stored.designoa$idVar)
putRcmdr("idpos", which(idlist %in% tclvalue(idVar)))
idEntry <- ttkcombobox(baseFrame, textvariable=idVar, width=20, values=idlist, state="readonly")
    tkbind(idEntry, "<<ComboboxSelected>>", idsel)
idHint <- ttklabel(baseFrame, text="(select, if desired)", foreground="#888888")

nrunVar <- tclVar(.stored.designoa$nrunVar)
putRcmdr("nrunpos", which(nrunslist %in% tclvalue(nrunVar)))
nrunEntry <- ttkcombobox(baseFrame, textvariable=nrunVar, width=5, values=nrunslist, state="readonly")
    tkbind(nrunEntry, "<<ComboboxSelected>>", nrunsel)
nrunHint <- ttklabel(baseFrame, text="(select, if desired)", foreground="#888888")

### allow for old version stored settings
if (is.null(.stored.designoa$optimVar)) optimVar <- tclVar("none")
else optimVar <- tclVar(.stored.designoa$optimVar)

putRcmdr("optimlist", c("none","min3","min34"))
putRcmdr("optimpos", which(optimlist %in% tclvalue(optimVar)))
optimEntry <- ttkcombobox(baseFrame, textvariable=optimVar, width=5, values=optimlist, state="readonly")

minrdfVar <- tclVar(.stored.designoa$minrdfVar)
minrdfEntry <- tkentry(baseFrame, width="8", textvariable=minrdfVar)
nrepVar <- tclVar(.stored.designoa$nrepVar)
nrepEntry <- tkentry(baseFrame, width="8", textvariable=nrepVar)
randomizeVariable <-  tclVar(.stored.designoa$cbInitials[2])
randomizecb <- ttkcheckbutton(baseFrame,text=gettextRcmdr("Randomization"),variable=randomizeVariable)
tkconfigure(randomizecb, takefocus=0)
seedVar <- tclVar(sample(31999,1))  ## always new
seedEntry <- tkentry(baseFrame, width="8", textvariable=seedVar)
tkconfigure(seedEntry, takefocus=0)
repeat.onlyVariable <- tclVar(.stored.designoa$cbInitials[1])
repeat.onlycb <- ttkcheckbutton(baseFrame,text=gettextRcmdr("Repeat only"),variable=repeat.onlyVariable)
tkconfigure(repeat.onlycb, takefocus=0)

inspectFrame <- ttklabelframe(tab1,text=gettextRcmdr("Inspect avaialable designs"))
fromVar <- tclVar(as.character(.stored.designoa$fromVar))  ## as.character should also make it work 
toVar <- tclVar(as.character(.stored.designoa$toVar))      ## with stored designs from previous version
fromRuns <- tkentry(inspectFrame, width="5", textvariable=fromVar)
toRuns <- tkentry(inspectFrame, width="5", textvariable=toVar)
tkgrid(tklabel(inspectFrame, text=gettextRcmdr("Number of of runs")), sticky="w", columnspan=4)
tkgrid(tklabel(inspectFrame, text="from"), fromRuns, tklabel(inspectFrame, text="to"), toRuns)

inspectButton <- buttonRcmdr(inspectFrame, text=gettextRcmdr("Show available designs\n     The menu remains open, \n     fetch it back after looking at designs"),
        foreground = "darkgreen", command = onInspect, 
        default = "normal", borderwidth = 3)
tkconfigure(inspectButton, takefocus=0)
parentcbVariable <- tclVar(.stored.designoa$cbInitials[9])
parentcb <- ttkcheckbutton(inspectFrame,variable=parentcbVariable,text=gettextRcmdr("Restrict to parent arrays"))
tkconfigure(parentcb, takefocus=0)
tkgrid(parentcb, columnspan=4, sticky="w")

usenlevelscbVariable <- tclVar(.stored.designoa$cbInitials[6])
usenlevelscb <- ttkcheckbutton(inspectFrame,variable=usenlevelscbVariable,text=gettextRcmdr("Use level information from factor details"))
tkconfigure(usenlevelscb, takefocus=0)
tkgrid(usenlevelscb, columnspan=4, sticky="w")

tkgrid(inspectButton, columnspan=4)

## preparations for bottom frame
bottomFrame <- tkframe(topdes2)

## grid base frame
## omitted nfaccb, on form, nfactors must always be specified
tkgrid(nfaclab <- tklabel(baseFrame, text=gettextRcmdr("Number of factors")), nfacEntry, nfacHint, sticky="w")
tkgrid(idlab <- tklabel(baseFrame, text=gettextRcmdr("Specific array")), idEntry, idHint, sticky="w")
tkgrid(optimlab <- tklabel(baseFrame, text=gettextRcmdr("Automatic optimization")), optimEntry, sticky="w")
tkgrid(nrunlab <- tklabel(baseFrame, text=gettextRcmdr("Minimum number of runs")), nrunEntry, nrunHint, sticky="w")

tkgrid(minrdflab <- tklabel(baseFrame, text=gettextRcmdr("Minimum residual df")), minrdfEntry, sticky="w", pady=c(15,0))
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
tkgrid.configure(nameFrame, pady=c(5,20))
tkgrid.configure(helptab1Button, sticky="ne")
tkgrid(baseFrame, inspectFrame, sticky="nw",columnspan=3)

## Factor Details Tab
## factor details frame
### facnameAutoVariable (not needed any more) and faclevelCommonVariable

putRcmdr("colnospecifyVariable", tclVar(.stored.designoa$cbInitials[3]))
colnoFrame <- tkframe(tab2)
colnospecifycb <- ttkcheckbutton(colnoFrame,text=gettextRcmdr("Manually specify column numbers for array ?"),
    variable=colnospecifyVariable, command=onRefresh)

tkgrid(colnospecifycb, sticky="w")
tkgrid(colnohint <- tklabel(colnoFrame, text=gettextRcmdr("useful for experts in customizing design properties\n(default column allocations are not always best)")),sticky="w")
if (idpos > 1) tkgrid(colnoFrame, pady=10,sticky="w")

## factor details
## values as vectors
    facnamlistt <- .stored.designoa$facnamlist
    nlevlistt <- .stored.designoa$nlevlist
    colnolistt <- .stored.designoa$colnolist
    faclevlistt <- .stored.designoa$faclevlist
    faclablistt <- .stored.designoa$faclablist
    varlistshortt <- if (as.numeric(tclvalue(nfacVar))<=50) 
                 Letters[1:tclvalue(nfacVar)] else paste("F",1:tclvalue(nfacVar),sep="")

    enterlistFrame <- ttkframe(tab2)
    listFrame <- ttklabelframe(enterlistFrame, text="Factor Details")
    putRcmdr("selpos", 1)
    putRcmdr("curfac", tclVar(varlistshortt[1]))
    putRcmdr("curfnam", tclVar(facnamlistt[1]))
    putRcmdr("curcolno", tclVar(colnolistt[1]))
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
    colno <- ttkentry(enterFrame, textvariable=curcolno, width=6)
    tkbind(colno, "<FocusOut>", colnochange)
    if (tclvalue(idVar)=="NULL") tkconfigure(colno, state="disabled")
    nlev <- ttkentry(enterFrame, textvariable=curnlev, width=6)
    tkbind(nlev, "<FocusOut>", nlevchange)
    flev <- ttkentry(enterFrame, textvariable=curflev, width=20)
    tkbind(flev, "<FocusOut>", flevchange)
    if (idpos>1 & tclvalue(colnospecifyVariable)=="1"){
    flab <- ttkentry(enterFrame, textvariable=curflab, width=15)
    } else
    flab <- ttkentry(enterFrame, textvariable=curflab, width=20)
    tkbind(flab, "<FocusOut>", flabchange)
    tkbind(flab, "<Key-Tab>", tabflab)
    if (idpos>1 & tclvalue(colnospecifyVariable)=="1"){
    tkgrid(tklabel(enterFrame,text=gettextRcmdr("Select"),width=6),
           tklabel(enterFrame,text=gettextRcmdr("Factor name"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("Column \nnumber"), width=6),
           tklabel(enterFrame,text=gettextRcmdr("no. of \nlevels"), width=6),
           tklabel(enterFrame,text=gettextRcmdr("Levels \n(separate with blanks)"), width=20),
           tklabel(enterFrame,text=gettextRcmdr("Comment/label \nfor html export"), width=15),
           sticky="w")
    tkgrid(fsel,fnam, colno, nlev, flev, flab, sticky="w")} else{
    tkgrid(tklabel(enterFrame,text=gettextRcmdr("Select"),width=6),
           tklabel(enterFrame,text=gettextRcmdr("Factor name"), width=15),
           tklabel(enterFrame,text=gettextRcmdr("no. of \nlevels"), width=6),
           tklabel(enterFrame,text=gettextRcmdr("Levels \n(separate with blanks)"), width=20),
           tklabel(enterFrame,text=gettextRcmdr("Comment or label \n(for html export only)"), width=20),
           sticky="w")
    tkgrid(fsel,fnam, nlev, flev, flab, sticky="w")
    }
    
    putRcmdr("facnamlist", tclVar(facnamlistt))
    putRcmdr("varlistshort", tclVar(varlistshortt))
    putRcmdr("nlevlist", tclVar(nlevlistt))
    putRcmdr("colnolist", tclVar(colnolistt))
    putRcmdr("faclevlist", tclVar(faclevlistt))
    putRcmdr("faclablist", tclVar(faclablistt))

    facshortListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=varlistshort,
        width = 6, background="#EBEBDC")
    tkbind(facshortListBox, "<<TraverseIn>>",function() tkfocus(fsel))

    facnameListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=facnamlist,
        width = 15, background="#EBEBDC")
    colnoListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=colnolist,
        width = 6, background="#EBEBDC")
    nlevListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=nlevlist,
        width = 6, background="#EBEBDC")
    faclevListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclevlist,
        width = 20, background="#EBEBDC")
    if (idpos>1 & tclvalue(colnospecifyVariable)=="1")
    faclabListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclablist,
        width = 15, background="#EBEBDC") else
    faclabListBox <- tklistbox(listFrame, height = min(10, as.numeric(tclvalue(nfacVar))), 
        selectmode = single, exportselection = "TRUE", listvariable=faclablist,
        width = 20, background="#EBEBDC")

    ## determine current index and reordering for onUp and onDown
    tkbind(facshortListBox, "<<ListboxSelect>>", checkIndexShort)
    tkbind(facnameListBox, "<<ListboxSelect>>", checkIndexNam)
    tkbind(nlevListBox, "<<ListboxSelect>>", checkIndexNlev)
    tkbind(colnoListBox, "<<ListboxSelect>>", checkIndexColno)
    tkbind(faclevListBox, "<<ListboxSelect>>", checkIndexFlev)
    tkbind(faclabListBox, "<<ListboxSelect>>", checkIndexLab)


    ### funktioniert, ist aber noch nicht schön
    scrollbar <- ttkscrollbar(listFrame, command = function(...) {
            tkyview(facshortListBox, ...)
            tkyview(facnameListBox, ...)
            tkyview(nlevListBox, ...)
            tkyview(colnoListBox, ...)
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

    if (idpos>1 & tclvalue(colnospecifyVariable)=="1")
    tkgrid(scrollbar, facshortListBox, facnameListBox, colnoListBox, nlevListBox, faclevListBox, faclabListBox, downupFrame, sticky = "nw") else
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
etyperbVariable <- tclVar(.stored.designoa$etyperbVariable)
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
decimalrbVariable <- tclVar(.stored.designoa$decimalrbVariable)
defaultrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("default"),variable=decimalrbVariable, value="default")
pointrb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr("."),variable=decimalrbVariable, value=".")
commarb <- tkradiobutton(decimalradioFrame,text=gettextRcmdr(","),variable=decimalrbVariable, value=",")
tkgrid(defaultrb, sticky="w")  ## in this case, leave default option from options
tkgrid(pointrb, sticky="w")
tkgrid(commarb, sticky="w")

## export directory
dirFrame <- ttklabelframe(tab6, text=gettextRcmdr("Storage Directory"))
putRcmdr("dirVar", tclVar(.stored.designoa$dirVar))
dirEntry <- tkentry(dirFrame, width="50", textvariable=dirVar)
dirButton <- buttonRcmdr(dirFrame, text = gettextRcmdr("Change directory"), 
        foreground = "darkgreen", width = "20", command = onChangeDir, 
        default = "normal", borderwidth = 3)
tkgrid(dirEntry, tklabel(dirFrame, text="   "), dirButton, sticky="w")

## export file name
putRcmdr("fileVar", tclVar(.stored.designoa$fileVar))
fileEntry <- tkentry(tab6, width="20", textvariable=fileVar)
efnamelabel <- tklabel(tab6,text=gettextRcmdr("Export file names: name below with appropriate endings (html or csv, and rda)"))
replacecbVariable <- tclVar(.stored.designoa$cbInitials[8])
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

OKCancelHelp(window=topdes2, helpSubject="Menu.oa")
tkconfigure(OKbutton, takefocus=0)
tkconfigure(cancelButton, takefocus=0)
tkconfigure(helpButton, takefocus=0)

tkgrid(buttonsFrame, bottomFrame, sticky="ew")


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
# End of Menu.oa
###
