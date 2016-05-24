aaggr0<-function () 
{
Library("YaleToolkit")
    defaults <- list(initial.variable = NULL, initial.mean = "0",initial.effectif="0")
    dialog.values <- getDialog("aaggr0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Donn", "\U00E9", "es manquantes",sep = "")))
    variableBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variable, 
            "all"))
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean","effectif"), initialValues = c(dialog.values$initial.mean,dialog.values$initial.effectif), 
        labels = gettextRcmdr(c("Rangement en fonction des effectifs","Effectifs (plutot que %)")))
  
 
    onOK <- function() {
        variables <- getSelection(variableBox)
        meanVar <- tclvalue(meanVariable)
        effectifVar <- tclvalue(effectifVariable)
        putDialog("aaggr0", list(initial.variable = variables, 
            initial.mean = meanVar,initial.effectif=effectifVar))
        closeDialog()
        if (length(variables) == 0) {
            errorCondition(recall = aaggr0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("plotMissing(", .activeDataSet, "[,c(\"", 
            listvar, "\")],tri=", meanVar, ",effectif=",effectifVar,")", sep = "")
        logger(command)
        justDoIt(command)

       command <- paste("statMissing(", .activeDataSet, "[,c(\"", 
            listvar, "\")],tri=", meanVar, ",effectif=",effectifVar,")", sep = "")
        doItAndPrint(command)

        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "is.na", reset = "aaggr0")
    tkgrid(getFrame(variableBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")

    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}