
bbigraphe0<-function(){
Library("qgraph")

    defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("bbigraphe0", defaults)
    initializeDialog(title = gettextRcmdr("Graphe des relations"))
    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "all"))
    onOK <- function() {
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 2) {
            errorCondition(recall = bbigraphe0, message = gettextRcmdr("You must select at least two variables"))
            return()
        }
        putDialog("bbigraphe0", list(initial.variables = variables))
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("bigraphe0(", .activeDataSet, "[c(\"", 
            listvar, "\")])", sep = "")
        logger(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot)
    OKCancelHelp(helpSubject = "qgraph", reset = "bbigraphe0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}

