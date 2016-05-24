Menu.EffectPlots <- function(){
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))
    putRcmdr("resp.list", response.names(eval(parse(text=.activeDataSet))))

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
      code <- as.logical(as.numeric(tclvalue(codeVar)))
      autolab <- as.logical(as.numeric(tclvalue(autolabVar)))
      alpha <- as.numeric(tclvalue(alphaVar))
      half <- as.logical(as.numeric(tclvalue(halfVar)))
      response <- getSelection(sel.resps)

     ## selected plots
       command <- paste("DanielPlot(", .activeDataSet, ", code=", as.character(code), ", autolab=", as.character(autolab), 
            ", alpha=", alpha, ", half=", as.character(half),", response=", dquote(response), ")", sep="")
       logger(command)
       hilf <- justDoItDoE(command)
       if (class(hilf)[1] == "try-error"){
             errorCondition(window=top,recall=NULL, message=gettextRcmdr(hilf))
             return()
       }
     if (length(grep("splitplot",di$type)) > 0)
         warning("Estimated effects for whole plot factors can stick out due to reasons not related to their importance!")
     closeDialog(window=top)
    }
    
    onalpha <- function(){
        oldwarn <- options("warn")$warn
        options(warn=0)
        alpha <- as.numeric(as.character(tclvalue(alphaVar)))
        options(warn=oldwarn)
        if (is.na(alpha)) 
           tk_messageBox(type="ok", message="Significance level must be numeric, please correct!", caption="Invalid alpha")
        else if (alpha<0 | alpha>1)
           tk_messageBox(type="ok", message="Significance level must be between 0 and 1, please correct!", caption="Invalid alpha")
    }
    onReset <- function(){
       putRcmdr("halfVar", tclVar("1"))
       putRcmdr("codeVar", tclVar("1"))
       putRcmdr("alphaVar", tclVar("0.1"))
       putRcmdr("autolabVar", tclVar("1"))
       tkconfigure(halfcb, variable=halfVar)
       tkconfigure(codecb, variable=codeVar)
       tkconfigure(autolabcb, variable=autolabVar)
       tkconfigure(alphaEntry, textvariable=alphaVar)
    }
   
   initializeDialog(title=gettextRcmdr("Effects plots for 2-level factors"))
   if (!exists("halfVar")) putRcmdr("halfVar", tclVar("1"))
   if (!exists("codeVar")) putRcmdr("codeVar", tclVar("1"))
   if (!exists("alphaVar")) putRcmdr("alphaVar", tclVar("0.1"))
   if (!exists("autolabVar")) putRcmdr("autolabVar", tclVar("1"))
   selFrame <- tkframe(top)
   putRcmdr("sel.resps", variableListBox(selFrame, variableList=resp.list, listHeight=10, 
        title="Response to be analysed (select one)",selectmode="single", initialSelection=0))
   tkgrid(selFrame,sticky="n")
   tkgrid(getFrame(sel.resps), sticky="n")

   resetBut <- tkbutton(top, text="Reset to defaults", command=onReset)
   halfcb <- ttkcheckbutton(top, text="Half normal plot ?", variable=halfVar)
   codecb <- ttkcheckbutton(top, text="Label effects with codes instead of names ?",
        variable=codeVar)
   alphaEntry <- ttkentry(top, textvariable=alphaVar,width="7")
   tkbind(alphaEntry,"<FocusOut>", onalpha)
   autolabcb <- ttkcheckbutton(top, text="Label significant effects only ?",
        variable=autolabVar)
   tkgrid(halfcb, resetBut, sticky="w")
   tkgrid.configure(resetBut, sticky="e")
   tkgrid(codecb, sticky="w")
   tkgrid(tklabel(top,text="Enter significance level for labelling:"), sticky="w")
   tkgrid(alphaEntry, sticky="w", padx="10")
   tkgrid(autolabcb, sticky="w")

    OKCancelHelp(helpSubject="Menu.EffectsPlots")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=2)
}