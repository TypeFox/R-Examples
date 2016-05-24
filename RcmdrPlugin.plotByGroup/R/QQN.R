QQN <- function () {
    defaults <- list(initial.group = NULL, initial.response = NULL)
    dialog.values <- getDialog("qqNormplots", defaults)
    initializeDialog(title = gettextRcmdr("qqNormplots by group"),
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
        putDialog ("qqNormplotByGroup", list (initial.group =groups,
                                          initial.response = response))
        closeDialog()
        if (length(groups) == 0) {
            errorCondition(recall = QQN,
                           message = gettextRcmdr(
                           "You must select at least one factor."))
            return()
        }
        if (length(response) == 0) {
            errorCondition(recall = QQN,
                           message = gettextRcmdr(
                           "You must select a response variable."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        groups.list <- paste(paste(groups, "=", .activeDataSet,
                                   "$", groups, sep =""),
                             collapse = ", ")
       panel <- ",prepanel = prepanel.qqmathline,
       panel = function(x, ...) {
          panel.qqmathline(x, ...)
          panel.qqmath(x, ...)
       }"
        doItAndPrint(paste("qqmath(~", response, " | ",
                           paste(groups, collapse = "*"), ", data=",
                           .activeDataSet,panel, ")",sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "qqmath", reset = "QQN")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 4, columns = 2)
}
