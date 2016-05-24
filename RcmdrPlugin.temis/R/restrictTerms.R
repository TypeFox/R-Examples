restrictTermsDlg <- function() {
    initializeDialog(title=.gettext("Select or Exclude Terms"))

    radioButtons(name="what",
                 buttons=c("retain", "exclude"),
                 labels=c(.gettext("Retain only these terms"),
                          .gettext("Exclude these terms")),
                 right.buttons=FALSE)

    tclTerms <- tclVar("")
    entryTerms <- ttkentry(top, width="30", textvariable=tclTerms)

    onOK <- function() {
        termsList <- strsplit(tclvalue(tclTerms), " ")[[1]]

        if(length(termsList) == 0) {
            .Message(.gettext("Please enter at least one term."), "error", parent=top)

            return()
        }
        else if(!all(termsList %in% colnames(dtm))) {
            wrongTerms <- termsList[!termsList %in% colnames(dtm)]
            .Message(sprintf(.ngettext(length(wrongTerms),
                                      "Term \'%s\' does not exist in the corpus.",
                                      "Terms \'%s\' do not exist in the corpus."),
                                      # TRANSLATORS: this should be opening quote, comma, closing quote
                            paste(wrongTerms, collapse=.gettext("\', \'"))),
                     "error", parent=top)

            return()
        }

        closeDialog()


        what <- tclvalue(whatVariable)
        if(what == "retain")
            doItAndPrint(paste("dtm <- dtm[, colnames(dtm) %in% c(\"",
                               paste(termsList, collapse="\", \""), "\")]", sep=""))
        else
            doItAndPrint(paste("dtm <- dtm[, !colnames(dtm) %in% c(\"",
                               paste(termsList, collapse="\", \""), "\")]", sep=""))

        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="restrictTermsDlg")
    tkgrid(whatFrame, sticky="w", columnspan=2, pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Terms (space-separated):")),
           columnspan=2, sticky="w")
    tkgrid(entryTerms, columnspan=2, sticky="w")
    tkgrid(buttonsFrame, sticky="ew", pady=6)
    dialogSuffix(focus=entryTerms)
}
