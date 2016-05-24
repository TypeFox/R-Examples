## one instance of assign replaced by justDoIt

Menu.param <- function(){
   designs <- listDesigns()
   
   
   onOK <- function(){
      newname <- tclvalue(nameVar)
      inner <- getSelection(innerBox)
      outr <- getSelection(outerBox)
      direction <- as.character(tclvalue(directionrbVar))

        if (is.element(newname, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(newname, gettextRcmdr("Object"))))
            {
              errorCondition(window=top,recall=Menu.param, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }
      
      command <- paste("param.design(inner=", inner, 
            ", outer=", outr, ", direction=\"",direction,"\")",sep="")
      hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=top,recall=Menu.param, message=gettextRcmdr(hilf))
             return()
            }
      logger(paste(newname, "<-", command))
        putRcmdr("hilf", hilf)
        ## replace assign by justDoIt; assign(newname, hilf, envir=.GlobalEnv)
        justDoIt(paste(newname, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
     closeDialog(window=top)
     activeDataSet(newname)
     tkfocus(CommanderWindow())
    }
    

    onname <- function(){
        if (!is.valid.name(tclvalue(nameVar))) 
           tk_messageBox(type="ok", message="invalid name for design, please correct!", caption="Invalid name")
    }   
   initializeDialog(title=gettextRcmdr("Taguchi parameter design ..."))
   nameVar <- tclVar("paramDesign.1")
   nameEntry <- ttkentry(top, textvariable=nameVar,width="25")
   
   innerBox <- variableListBox(top, variableList=designs, selectmode="single",
        title=gettextRcmdr("Inner array (select one design)"),
        initialSelection=NULL)
   outerBox <- variableListBox(top, variableList=designs, selectmode="single",
        title=gettextRcmdr("Outer array (select one design)"),
        initialSelection=NULL)

   directionrbVar <- tclVar("long")
   longrb <- ttkradiobutton(top, text="Long version ",
        variable=directionrbVar,value="long")
   widerb <- ttkradiobutton(top, text="Wide version (crossed, like usual in Taguchi representation)",
        variable=directionrbVar,value="wide")
   
   tkbind(nameEntry,"<FocusOut>", onname)
   tkgrid(tklabel(top,text="Enter name for new design:"), sticky="w")
   tkgrid(nameEntry, sticky="w",padx="20")
   tkgrid(tklabel(top, text="Select two different designs with non-overlapping factor names:"), sticky="w")
   tkgrid(getFrame(innerBox), getFrame(outerBox), sticky="n")
   tkgrid(longrb, sticky="w")
   tkgrid(widerb, sticky="w")

    OKCancelHelp(helpSubject="Menu.param")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=2)
}