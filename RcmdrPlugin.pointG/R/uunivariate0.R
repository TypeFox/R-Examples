uunivariate0 <- function (){
    defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("uunivariate0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Analyse graphique univari", "\U00E9", "e",sep = "")))

    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variables, 
            "all"))
    onOK <- function() {
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 1) {
            errorCondition(recall = uunivariate0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        putDialog("uunivariate0", list(initial.variables = variables))
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("univariate0(", .activeDataSet, "[c(\"", 
            listvar, "\")])", sep = "")
        logger(command)
        justDoIt(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot)
    OKCancelHelp(helpSubject = "univariate0", reset = "uunivariate0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}
