Menu.cand <- function(){
    ### menu for selecting or creating candidate design for D-optimization
    initializeDialog(window=top, title=gettextRcmdr("Which type of candidate design ?"))
    

    changeActive <- function(){
         closeDialog(window=top)
         Menu.import()
    }

    changeActiveAll <- function(){
         closeDialog(window=top)
         selectActiveDataSet()
    }
    selectActiveDataSet <- function(){
       ## copied from Rcmdr data menu
       dataSets <- listDataSets()
      .activeDataSet <- ActiveDataSet()
      if ((length(dataSets) == 1) && !is.null(.activeDataSet)) {
        Message(message=gettextRcmdr("There is only one dataset in memory."),
         type="warning")
        tkfocus(CommanderWindow())
        return()
      }
      if (length(dataSets) == 0){
        Message(message=gettextRcmdr("There are no data sets from which to choose."),
          type="error")
        tkfocus(CommanderWindow())
        return()
      }
      initializeDialog(title=gettextRcmdr("Select Data Set"))
      dataSetsBox <- variableListBox(top, dataSets, title=gettextRcmdr("Data Sets (pick one)"),
        initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
      onOK <- function(){
        activeDataSet(getSelection(dataSetsBox))
        closeDialog()
        tkfocus(CommanderWindow())
      }
      OKCancelHelp()
      tkgrid(getFrame(dataSetsBox), sticky="nw")
      tkgrid(buttonsFrame, sticky="w")
      dialogSuffix(rows=2, columns=1)
    }

    newff <- function(){
         closeDialog(window=top)
         Menu.fac()
    }
    newoa <- function(){
         closeDialog(window=top)
         Menu.oa()
    }
    newlhs <- function(){
         closeDialog(window=top)
         Menu.lhs()
    }
    newFrF2 <- function(){
         closeDialog(window=top)
         Menu.FrF2level()
    }
    onCancel <- function() {
        if (GrabFocus()) 
            tkgrab.release(top)
        tkdestroy(top)
    }
    onOK <- onCancel
    
    if (is.null(ActiveDataSet())) tkgrid(tklabel(top, text=gettextRcmdr("Per default, the active data frame is used as the candidate design.\nCurrent active data frame: NONE")), sticky="w")
    else {tkgrid(tklabel(top, text=gettextRcmdr("Candidate design (active data frame):")),sticky="w")
          tkgrid(tklabel(top, text=ActiveDataSet()),sticky="w")}
    tkgrid(tklabel(top, text=gettextRcmdr("Select new candidate data frame (will become active data frame)")),sticky="w")
    tkgrid(tkbutton(top, text=gettextRcmdr("from existing designs"), command=changeActive, width="100"), sticky="we", pady="5")
    tkgrid(tkbutton(top, text=gettextRcmdr("from existing data frames in general"), command=changeActiveAll, width="100"), sticky="we", pady="5")
    tkgrid(ttkseparator(top, orient="horizontal"),pady="10")
    tkgrid(tklabel(top, text=gettextRcmdr("Create new candidate design\nThe buttons below are quick links to the most frequently used types of candidate designs."),justify="left"),sticky="w",pady="5")
    tkgrid(tkbutton(top, text=gettextRcmdr("Full factorial design"), command=newff, width="100"), sticky="we", pady="5")
    tkgrid(tkbutton(top, text=gettextRcmdr("Orthogonal array at arbitrary numbers of levels"), command=newoa, width="100"), sticky="we", pady="5")
    tkgrid(tkbutton(top, text=gettextRcmdr("Latin hypercube design with only quantitative factors"), command=newlhs, width="100"), sticky="we", pady="5")
    tkgrid(tkbutton(top, text=gettextRcmdr("Fractional factorial 2-level design"), command=newFrF2, width="100"), sticky="we", pady="5")
    
        cancelButton <- buttonRcmdr(top, text = gettextRcmdr("Cancel"), 
            foreground = "red", width = "12", command = onCancel, 
            borderwidth = 3)
    tkgrid(cancelButton, sticky="s",pady="20")
    dialogSuffix(window=top, rows=5, columns=1)
}
