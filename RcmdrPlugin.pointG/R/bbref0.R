bbref0<-function () 
{
    defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("bbref0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Bref r", "\U00E9", "sum","\U00E9",sep = "")))
    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variables, 
            "all"))
    onOK <- function() {
        variables <- getSelection(variablesBox)
        closeDialog()
        if (length(variables) < 1) {
            errorCondition(recall = bbref0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        putDialog("bbref0", list(initial.variables = variables))
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("bref.df0(", .activeDataSet, "[c(\"", 
            listvar, "\")])", sep = "")
        logger(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    groupsBox(scatterPlot)
    OKCancelHelp(helpSubject = "bref.variable0", reset = "bbref0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}