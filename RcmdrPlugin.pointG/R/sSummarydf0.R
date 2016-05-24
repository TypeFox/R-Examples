sSummarydf0<-
function () 
{
    defaults <- list(initial.variables = NULL,initial.effectif="1")
    dialog.values <- getDialog("sSummarydf0", defaults)
    initializeDialog(title = gettextRcmdr("Resume numerique univarie"))
    variablesBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "all"))
 checkBoxes(frame = "checkBoxFrame", boxes = c("effectif"), initialValues = c(dialog.values$initial.effectif), 
        labels = gettextRcmdr(c("% (plutot qu'en effectif) pour les variables categorielles")))

    onOK <- function() {
        variables <- getSelection(variablesBox)
        effectifVar <- tclvalue(effectifVariable)
     putDialog("sSummarydf0", list(initial.variables = variables,initial.effectif=effectifVar))
   
        closeDialog()
        if (length(variables) < 1) {
            errorCondition(recall = sSummarydf0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("Summarydf0(", .activeDataSet, "[c(\"", 
            listvar, "\")],pourcent=",effectifVar,")", sep = "")
        logger(command)
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    
    OKCancelHelp(helpSubject = "Summaryvariable0", reset = "sSummarydf0")
    tkgrid(getFrame(variablesBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")
   tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 6, columns = 2)
}