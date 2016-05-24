ccroise0<-function () 
{
Library("gpairs")
    defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("ccroise0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Graphiques crois", "\U00E9", "s", 
            sep = "")))


    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variables, 
            "all"))
    onOK <- function() {
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 2) {
            errorCondition(recall = ccroise0, message = gettextRcmdr("You must select at least two variables"))
            return()
        }
        putDialog("ccroise0", list(initial.variables = variables))
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("gpairs(", .activeDataSet, "[c(\"", 
            listvar, "\")])", sep = "")
        logger(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot)
    OKCancelHelp(helpSubject = "gpairs", reset = "ccroise0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}
