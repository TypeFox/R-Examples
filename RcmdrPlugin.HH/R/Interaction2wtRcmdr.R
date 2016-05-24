"Interaction2wtRcmdr" <-
function() {
    if (length(grep("HH", search()))==0) stop("Please attach the HH directory.")
    initializeDialog(title=gettextRcmdr("Interaction twoway table"))
    groupBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factors (pick two or more)"), selectmode="multiple")
    responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Response Variable (pick one)"))
    onOK <- function() {
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
        closeDialog()
        if (2 > length(groups)) {
            errorCondition(recall=Interaction2wtRcmdr, message=gettextRcmdr("Select at least two factors."))
            return()
            }
        if (0 == length(response)) {
            errorCondition(recall=Interaction2wtRcmdr, message=gettextRcmdr("No response variable selected."))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        i2wt.command <- paste("interaction2wt(",
                              response, " ~ ", paste(groups, collapse=' + '),
                              ", data=", .activeDataSet,
                              if ("1" == tclvalue(simpleVariable)) ", simple=TRUE",
                              ')', sep="")
        doItAndPrint(i2wt.command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    optionsFrame <- tkframe(top)
    buttonsFrame <- tkframe(top)
    OKCancelHelp(helpSubject="interaction2wt")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky="nw")

    simpleVariable <- tclVar("0")
    simpleFrame <- tkframe(top)
    simpleCheckBox <- tkcheckbutton(simpleFrame, variable=simpleVariable)
    tkgrid(tklabel(simpleFrame, text=gettextRcmdr("Simple Effects")),
           simpleCheckBox, sticky="w")
    tkgrid(simpleFrame, sticky="w")

    tkgrid(optionsFrame, columnspan=2, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
  }

