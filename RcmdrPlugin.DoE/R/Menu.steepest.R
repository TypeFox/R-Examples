Menu.steepest <- function(){
     ## menu for steepest slope (ascent or descent, starting at 0,0,0,0 or 
     ## at the stationary point
     
     ## dist should also be on offer
     ## currently uses the default only
     ## default and acceptable values for dist should depend on type
    .activeModel <- ActiveModel()
    degree <- get(.activeModel)$order
    initializeDialog(window=top, title=gettextRcmdr("Steepest slope analysis"))
    
    if (degree==2) typerbVariable <- tclVar("canonical.path") 
        else typerbVariable <- tclVar("steepest")
    
    if (tclvalue(typerbVariable)=="steepest") 
        explainVariable <- tclVar("positive values only, separate by blank\nsteps from all xs at 0")
    else
        explainVariable <- tclVar("positive values only, separate by blank\none direction from all xs at 0 or\n+/- directions from stationary point")
    if (degree==2){
       typerbFrame <- ttklabelframe(top,text=gettextRcmdr("Where to start steepest slope ?"))
       midrbButton <- tkradiobutton(typerbFrame, text="steepest slope from all xs at 0",variable=typerbVariable, value="steepest")
       statrbButton <- tkradiobutton(typerbFrame, text="steepest slope from stationary point",variable=typerbVariable, value="canonical.path")
       tkgrid(midrbButton, sticky ="w")
       tkgrid(statrbButton, sticky ="w")
    }
    
    dirrbVariable <- tclVar("ascent")
    dirrbFrame <- ttklabelframe(top,text=gettextRcmdr("Up or Down ?"))
    uprbButton <- tkradiobutton(dirrbFrame, text="steepest ascent",variable=dirrbVariable, value="ascent")
    downrbButton <- tkradiobutton(dirrbFrame, text="steepest descent",variable=dirrbVariable, value="descent")
    tkgrid(uprbButton, sticky ="w")
    tkgrid(downrbButton, sticky ="w")
    if (degree==2) tkgrid(typerbFrame, dirrbFrame, sticky="w")
    else tkgrid(ttklabel(top, text="Steepest slope starts at all xs at 0"), dirrbFrame, sticky="w")
    
    if (!exists(".move.distances", where="RcmdrEnv")) 
         putRcmdr(".move.distances", seq(0,5,0.5))
    distVar <- tclVar(paste(.move.distances, collapse=" "))
    distEntry <- tkentry(top, width="50", textvariable=distVar)
    distlab <- ttklabel(top, text="Step widths for steepest directions")
    #if (!degree==2 || tclvalue(typerbVariable)=="steepest")
    distexplain <- ttklabel(top, textvariable=explainVariable)
    #else
    #distexplain <- ttklabel(top, text="positive values only, separate by blank\n+/- directions from stationary point")
    tkgrid(distlab, distEntry,sticky="w")
    tkgrid(tklabel(top, text="NOTE:"), distexplain,sticky="e")
    tkgrid.configure(distlab, sticky="e")
    tkgrid.configure(distexplain, sticky="w")

    onOK <- function(){
        dist <- as.numeric(strsplit(tclvalue(distVar), " ", fixed=TRUE)[[1]])
        if (tclvalue(typerbVariable)=="canonical.path") 
           dist <- sort(unique(c(-dist, dist)))
        command <- paste(tclvalue(typerbVariable), "(",.activeModel,", descent=", tclvalue(dirrbVariable)=="descent", 
            ", dist=c(", paste(dist, collapse=","), "))")
        hilf <- doItAndPrint(command)
        if (class(hilf)[1]=="try-error") {
            errorCondition(window=top,recall=Menu.steepest, message=gettextRcmdr(hilf))
             return()
            }
        closeDialog(window=top)
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.steepest")
    tkgrid(buttonsFrame, sticky="w", columnspan=2)
    dialogSuffix(window=top, rows=2, columns=2)
    }
