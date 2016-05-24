Menu.loadcatlg <- function(){
    ## FrF2.catlg128 is loaded
    candlist <- ""
    logger("require(FrF2.catlg128)")
    if (require(FrF2.catlg128)){
         candlist <- c("catlg128.8to15",paste("catlg128",16:23,sep="."))}
    if (identical(candlist, "")) {
        Message(message=gettextRcmdr("There are no catalogues from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    
    initializeDialog(title=gettextRcmdr("Select catalogue of regular fractional factorial 2-level designs"))
    catlgBox <- variableListBox(top, candlist, title=gettextRcmdr("Catalogues (pick one)"),
        initialSelection=NULL)
## create button that calls further menu
## cascade menu in the menu entry
    onOK <- function(){
        command <- paste("data(",getSelection(catlgBox),")")
        logger(command)
        justDoItDoE(command)
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.loadcatlg")
    tkgrid(getFrame(catlgBox), sticky="n")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }
