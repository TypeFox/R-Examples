Menu.import <- function(){
    designs <- listDesigns()
    .activeDataSet <- ActiveDataSet()
  #  if ((length(designs) == 1) && !is.null(.activeDataSet) && designs==.activeDataSet) {
  #      Message(message=gettextRcmdr("There is only one design in memory."),
  #              type="warning")
  #      tkfocus(CommanderWindow())
  #      return()
  #      }
    if (length(designs) == 0){
        Message(message=gettextRcmdr("There are no designs from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    initializeDialog(title=gettextRcmdr("Select Design"))
    designsBox <- variableListBox(top, designs, title=gettextRcmdr("Designs (pick one)"),
        initialSelection=if (is.null(.activeDataSet) || !.activeDataSet %in% designs) NULL 
             else which(.activeDataSet == designs) - 1)
    onOK <- function(){
        activeDataSet(getSelection(designsBox))
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.import")
    tkgrid(getFrame(designsBox), sticky="n")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }
