barCode2 <- function () {
    Library("barcode")
defaults <- list(initial.x = NULL,initialGroup=NULL) 
    dialog.values <- getDialog("barCode2", defaults)
    initializeDialog(title = gettextRcmdr("Graphe en code barre"))
    xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
                            initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    

    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    onOK <- function() {
        x <- getSelection(xBox)
        
        putDialog ("barCode2", list(initial.x = x, 
                                   initial.group=if (.groups == FALSE) NULL else .groups))
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = barCode2, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        varx <- paste(.activeDataSet, "$", x, sep = "")
       
    if (is.null(.groups) || .groups == FALSE) {
            command <- (paste("barcode(", varx, ",xlab=\"",x,"\")", sep = ""))
            doItAndPrint(command)
        }
        else {
varg <- paste(.activeDataSet, "$", .groups, sep = "")

            command <- (paste("barcode(split(",varx,",",varg,"),xlab=\"",x,"\")", sep = ""))
            doItAndPrint(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(boxPlot, initialGroup=initial.group, 
              initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
    OKCancelHelp(helpSubject = "barcode", reset = "barCode2")
    tkgrid(getFrame(xBox), sticky = "nw")
    
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 4, columns = 1)
}

