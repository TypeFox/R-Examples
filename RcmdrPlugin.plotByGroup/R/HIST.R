HIST <- function () {
    defaults <- list(initial.group = NULL, initial.response = NULL)
    dialog.values <- getDialog("Histograms", defaults)
    initializeDialog(title = gettextRcmdr("Histograms by group"),
                     window=top)
    groupBox <- variableListBox(top, Factors(), selectmode = "multiple",
                                title = gettextRcmdr("Factors (pick one or more)"),
                                initialSelection =
                                varPosn(dialog.values$initial.group,
                                        "factor"))
    responseBox <- variableListBox(top,
                                   Numeric(),
                                   title = gettextRcmdr(
                                   "Response Variable(pick one)"),
                                   initialSelection =
                                   varPosn(dialog.values$initial.response,
                                           "numeric"))
    onOK <- function() {
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
        putDialog ("HistogramByGroup", list (initial.group =groups,
                                          initial.response = response))
        closeDialog()
        if (length(groups) == 0) {
            errorCondition(recall = HIST,
                           message = gettextRcmdr(
                           "You must select at least one factor."))
            return()
        }
        if (length(response) == 0) {
            errorCondition(recall = HIST,
                           message = gettextRcmdr(
                           "You must select a response variable."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        groups.list <- paste(paste(groups, "=", .activeDataSet,
                                   "$", groups, sep =""),
                             collapse = ", ")
        doItAndPrint(paste("histogram(~", response, " | ",
                           paste(groups, collapse = "*"), ", data=",
                           .activeDataSet, ")", sep = ""))
        #doItAndPrint("summary(fit)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "histogram", reset = "HIST")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 4, columns = 2)
}
