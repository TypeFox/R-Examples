## one instance of assign replaced by justDoIt

Menu.addcenter <- function(){
   .activeDataSet <- ActiveDataSet()
   di <- design.info(eval(parse(text=.activeDataSet)))
     if (length(grep("splitplot",di$type)) > 0){
         tk_messageBox(type="ok", message="Center points are not implemented for split plot designs.", caption="This does not work")
         return()
         }
   
   onOK <- function(){
      newname <- tclvalue(nameVar)
      ncenter <- as.numeric(tclvalue(ncenterVar))
      distribute <- as.numeric(tclvalue(distributeVar))

        if (is.element(newname, listObjects()))
          {
          if ("no" == tclvalue(checkReplace(newname, gettextRcmdr("Object"))))
            {
              errorCondition(window=top,recall=Menu.addcenter, 
              gettextRcmdr("Introduce another name for the new data.frame, or allow replacing."))
              return()
             }
          }

      command <- paste("add.center(", .activeDataSet, 
            ", ncenter=", ncenter, ", distribute=", distribute,")")
      
      hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=top,recall=Menu.addcenter, message=gettextRcmdr(hilf))
             return()
            }
      logger(paste(newname, "<-", command))
      ## replace assign by justDoIt; assign(newname, hilf, envir=.GlobalEnv)
      putRcmdr("hilf", hilf)
      justDoIt(paste(newname, "<- getRcmdr(\"hilf\")"))
        rm("hilf", pos="RcmdrEnv")
     closeDialog(window=top)
     activeDataSet(newname)
     tkfocus(CommanderWindow())
    }
    
    onncenter <- function(){
        oldwarn <- options("warn")$warn
        options(warn=0)
        ncenter <- as.numeric(as.character(tclvalue(ncenterVar)))
        options(warn=oldwarn)
        if (is.na(ncenter)) 
           tk_messageBox(type="ok", message="number of center points must be numeric, please correct!", caption="Invalid ncenter")
        else{ 
        if (ncenter<1)
           tk_messageBox(type="ok", message="number of center points must be between at least 1, please correct!", caption="Invalid ncenter")
        if (!ncenter%%1==0)
           tk_messageBox(type="ok", message="number of center points must be integer, please correct!", caption="Invalid ncenter")}
    }


    ondistribute <- function(){
        oldwarn <- options("warn")$warn
        options(warn=0)
        distribute <- as.numeric(as.character(tclvalue(distributeVar)))
        options(warn=oldwarn)
        if (is.na(distribute)) 
           tk_messageBox(type="ok", message="distribute must be numeric, please correct!", caption="Invalid distribute")
        else {
        if (distribute<1 | distribute>as.numeric(tclvalue(ncenterVar)))
           tk_messageBox(type="ok", message="distribute must be between 1 and ncenter, please correct!", caption="Invalid distribute")
        if (!distribute%%1==0)
           tk_messageBox(type="ok", message="distribute must be integer, please correct!", caption="Invalid distribute")
        }
    }
    onname <- function(){
        if (!is.valid.name(tclvalue(nameVar))) 
           tk_messageBox(type="ok", message="invalid name for design, please correct!", caption="Invalid name")
    }   
   initializeDialog(title=gettextRcmdr("Add center points to 2-level design"))
   nameVar <- tclVar(paste(.activeDataSet, "withcenterpts", sep="."))
   ncenterVar <- tclVar("6")
   distributeVar <- tclVar("3")
   nameEntry <- ttkentry(top, textvariable=nameVar,width="25")
   ncenterEntry <- ttkentry(top, textvariable=ncenterVar,width="5")
   distributeEntry <- ttkentry(top, textvariable=distributeVar,width="5")
   tkbind(nameEntry,"<FocusOut>", onname)
   tkbind(ncenterEntry,"<FocusOut>", onncenter)
   tkbind(distributeEntry,"<FocusOut>", ondistribute)
   tkgrid(tklabel(top,text="Enter name for new design:"), sticky="w")
   tkgrid(nameEntry, sticky="w", padx="20")
   tkgrid(tklabel(top,text="Enter number of center points:"), sticky="w")
   tkgrid(ncenterEntry, sticky="w", padx="20")
   tkgrid(tklabel(top,text="Enter number of positions over which center points are systematically distributed:"), sticky="w")
   tkgrid(distributeEntry, sticky="w", padx="20")

    OKCancelHelp(helpSubject="Menu.addcenter")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
}