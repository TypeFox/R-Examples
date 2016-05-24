nettoyage0<-
function () 
{
Library("YaleToolkit")
    defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("nettoyage0", defaults)
    initializeDialog(title = gettextRcmdr("Que sais-je"))
    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "all"))
    onOK <- function() {
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 1) {
            errorCondition(recall = nettoyage0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        putDialog("nettoyage0", list(initial.variables = variables))
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("whatis(", .activeDataSet, "[c(\"", 
            listvar, "\")])", sep = "")
        logger(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot)
    OKCancelHelp(helpSubject = "whatis", reset = "nettoyage0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}