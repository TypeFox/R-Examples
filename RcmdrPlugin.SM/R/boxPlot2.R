boxPlot2 <- function () {
    defaults <- list(initial.x = NULL, initial.identify = "none", initialGroup=NULL) 
    dialog.values <- getDialog("boxPlot2", defaults)
    initializeDialog(title = gettextRcmdr("Boxplot"))
    xBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variable (pick one)"), 
                            initialSelection = varPosn (dialog.values$initial.x, "numeric"))
    radioButtons(name = "identify", buttons = c("y", "identify", "none"), 
                 labels = gettextRcmdr(c("Automatique", "Avec la souris", "Non")), 
                 title = gettextRcmdr("Identification des points extr\xEAmes"), 
                 initialValue = dialog.values$initial.identify)
    initial.group <- dialog.values$initial.group
    .groups <- if (is.null(initial.group)) FALSE else initial.group
    onOK <- function() {
        x <- getSelection(xBox)
        identifyPoints <- tclvalue(identifyVariable)
        putDialog ("boxPlot", list(initial.x = x, initial.identify = identifyPoints, 
                                   initial.group=if (.groups == FALSE) NULL else .groups))
        closeDialog()
        if (length(x) == 0) {
            errorCondition(recall = boxPlot2, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        var <- paste(.activeDataSet, "$", x, sep = "")
        if (identifyPoints == "identify")
            RcmdrTkmessageBox(title = "Identify Points", 
                              message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
                                              gettextRcmdr(if (MacOSXP()) "esc key to exit."
                                                           else "right button to exit."), sep = ""), 
                              icon = "info", type = "ok")
        if (is.null(.groups) || .groups == FALSE) {
            command <- paste("Boxplot( ~ ", x, ", col=\"pink\",data=", .activeDataSet, ', id.method="', 
                             identifyPoints, '")', sep="")
            doItAndPrint(command)
        }
        else {
            command <- paste("Boxplot(", x, "~", .groups, ",col=\"pink\",data=", .activeDataSet, 
                             ', id.method="', identifyPoints, '")', sep = "")
            doItAndPrint(command)
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(boxPlot, initialGroup=initial.group, 
              initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
    OKCancelHelp(helpSubject = "boxplot", reset = "boxPlot2")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(identifyFrame, stick = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 4, columns = 1)
}

