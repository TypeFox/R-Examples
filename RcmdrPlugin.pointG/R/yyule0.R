yyule0<-function () 
{
    defaults <- list(initial.group = NULL, initial.variable = NULL, 
        initial.mean = "0")
    dialog.values <- getDialog("yyule0", defaults)
    initializeDialog(title = gettextRcmdr(paste("Profils de modalit", 
        "\U00E9", "s (Q)", sep = "")))
    groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Variable explicative categorielle"), 
        initialSelection = varPosn(dialog.values$initial.group, 
            "factor"))
    variableBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variables (select one or more)"), 
        selectmode = "multiple", initialSelection = varPosn(dialog.values$initial.variable, 
            "all"))
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean"), initialValues = c(dialog.values$initial.mean), 
        labels = gettextRcmdr(c("Rangement en fonction de la valeur de Q")))
    onOK <- function() {
        group <- getSelection(groupBox)
        variables <- getSelection(variableBox)
        meanVar <- tclvalue(meanVariable)
        putDialog("yyule0", list(initial.group = group, initial.variable = variables, 
            initial.mean = meanVar))
        closeDialog()
        if (length(variables) <2) {
            errorCondition(recall = yyule0, message = gettextRcmdr("You must select at least two variables"))
            return()
        }
        if (length(group) == 0) {
            errorCondition(recall = yyule0, message = gettextRcmdr("You must select a variable"))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        listvar <- paste(variables, collapse = "\",\"")
        names.levels <- eval(parse(text = paste("levels(", .activeDataSet, 
            "$", group, ")", sep = "")), envir = .GlobalEnv)
        nvalues <- length(names.levels)
        for (i in 1:nvalues) {
            namelevel <- names.levels[i]
            command <- paste("Yule04(", .activeDataSet, "$", 
                group, ",", .activeDataSet, "[,c(\"", listvar, 
                "\")],\"", namelevel, "\",nameYY=c(\"", listvar, 
                "\")\n", ",tri=", meanVar, ")", sep = "")
            doItAndPrint(command)
            command <- paste("Yule03(", .activeDataSet, "$", 
                group, ",", .activeDataSet, "[,c(\"", listvar, 
                "\")],\"", namelevel, "\",nameYY=c(\"", listvar, 
                "\")\n", ",tri=", meanVar, ")", sep = "")
            doItAndPrint(command)
            if (i < nvalues) {
                justDoIt("dev.new()")
            }
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "Yule03", reset = "yyule0")
    tkgrid(getFrame(groupBox), getFrame(variableBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 2, columns = 1)
}
