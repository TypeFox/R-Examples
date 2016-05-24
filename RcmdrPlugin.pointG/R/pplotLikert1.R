pplotLikert1<-function () 
{
Library("lattice")

    defaults <- list(initial.group = NULL, initial.variable = NULL, 
        initial.mean = "0",initial.scale="0")
    dialog.values <- getDialog("pplotLikert1", defaults)
    initializeDialog(title = gettextRcmdr("Echelles de Likert en saucisson"))
    groupBox <- variableListBox(top, Factors(), title = gettextRcmdr(paste("Variable explicative cat", "\U00E9", "gorielle",sep = "")), 
        initialSelection = varPosn(dialog.values$initial.group, 
            "factor"))
    variableBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variable, 
            "numeric"))
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean","scale"), initialValues = c(dialog.values$initial.mean,dialog.values$initial.scale), 
        labels = gettextRcmdr(c("Rangement en fonction de la moyenne (graphe) et de la p-value (statistiques)",paste("Mise ","\U00E0"," l'", "\U00E9", "chelle des abscisses",sep = ""))))
    onOK <- function() {
        group <- getSelection(groupBox)
        variables <- getSelection(variableBox)
        meanVar <- tclvalue(meanVariable)
        scaleVar <- tclvalue(scaleVariable)

        putDialog("pplotLikert1", list(initial.group = group, 
            initial.variable = variables, initial.mean = meanVar,initial.scale=scaleVar))
        closeDialog()
        if (length(variables) <2) {
            errorCondition(recall = pplotLikert1, message = gettextRcmdr("You must select at least two Likert variables"))
            return()
        }
        if (length(group) == 0) {
            errorCondition(recall = pplotLikert1, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        command <- paste("plotLikert1(", .activeDataSet, "[,c(\"", 
            listvar, "\")],", .activeDataSet, "$", group, ",tri=", 
            meanVar,",scale=",scaleVar, ")", sep = "")
        doItAndPrint(command)
        command <- paste("statLikert1(", .activeDataSet, "[,c(\"", 
            listvar, "\")],", .activeDataSet, "$", group, ",tri=", 
            meanVar, ")", sep = "")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "plotLikert1", reset = "pplotLikert1")
    tkgrid(getFrame(groupBox), getFrame(variableBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}
