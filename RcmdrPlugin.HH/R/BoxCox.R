"BoxCox" <-
function(){
    initializeDialog(title="Box-Cox Transformations")
    variablesBox <- variableListBox(top, Numeric(), selectmode="multiple",
        title="Select variables (one or more)")
    onOK <- function(){
        variables <- getSelection(variablesBox)
        if (length(variables) < 1) {
            errorCondition(recall=BoxCox,
              message="You must select one or more variables.")
            return()
            }
        closeDialog()
        command <- paste("box.cox.powers(na.omit(cbind(",
            paste(paste(variables, "=", ActiveDataSet(), "$", variables, sep=""),
                collapse=", "), ")))", sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="box.cox.powers")
    tkgrid(getFrame(variablesBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }

