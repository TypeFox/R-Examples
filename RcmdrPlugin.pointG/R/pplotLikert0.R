pplotLikert0<-
function () 
{
    defaults <- list(initial.variable = NULL, initial.mean = "0",initial.adaptation="0")
    dialog.values <- getDialog("pplotLikert0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Collection d'", "\U00E9", "chelles de Likert", 
            sep = "")))
    variableBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variable, 
            "numeric"))
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean","adaptation"), 
initialValues = c(dialog.values$initial.mean,dialog.values$initial.adaptation), 
        labels = gettextRcmdr(c("Rangement en fonction de la moyenne",paste("Adapter l'", "\U00E9", "chelle des abscisses", 
            sep = ""))))
    onOK <- function() {
        variables <- getSelection(variableBox)
        meanVar <- tclvalue(meanVariable)
	adaptationVar <- tclvalue(adaptationVariable)
        putDialog("pplotLikert0", list(initial.variable = variables, 
            initial.mean = meanVar,initial.adaptation=adaptationVar))
        closeDialog()
        if (length(variables) < 2) {
            errorCondition(recall = aaggr0, message = gettextRcmdr("Il faut au moins selectionner deux variables"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("plotLikert0(", .activeDataSet, "[,c(\"", 
            listvar, "\")],tri=", meanVar, ",adaptation=",adaptationVar,")", sep = "")
        logger(command)
        justDoIt(command)
        command <- paste("statLikert0(", .activeDataSet, "[,c(\"", 
            listvar, "\")],tri=", meanVar, ")", sep = "")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "plotLikert0", reset = "pplotLikert0")
    tkgrid(getFrame(variableBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}
