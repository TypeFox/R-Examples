termChisqDist <- function(term, dtm, n=5, variable=NULL) {
    if(!term %in% colnames(dtm))
         stop("'term' is not present in 'dtm'")

    dev <- sweep(as.matrix(dtm)/col_sums(dtm), 1,
                 as.matrix(dtm[, term])/sum(dtm[, term]), "-")
    chisq <- sweep(dev^2, 1, row_sums(dtm)/sum(dtm), "/")

    # na.rm=TRUE is here because some empty documents might exist, even if it's useless
    if(is.null(variable)) {
        head(sort(colSums(chisq, na.rm=TRUE)), n)
    }
    else {
        sapply(levels(factor(variable)), function(l) {
            if(sum(dtm[variable == l, term]) == 0)
                NA
            else
                head(sort(colSums(chisq[variable == l, , drop=FALSE], na.rm=TRUE)), n)
        })
    }
}

termCoocDlg <- function() {
    if(!(exists("dtm") && class(dtm) == "DocumentTermMatrix")) {
        .Message(message=.gettext("Please import a corpus and create the document-term matrix first."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Terms Co-occurring With Chosen Terms"))

    tclTerms <- tclVar("")
    entryTerms <- ttkentry(top,  width="35", textvariable=tclTerms)

    tclN <- tclVar(10)
    spinN <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                      inc=1, textvariable=tclN,
                      validate="all", validatecommand=.validate.uint)

    vars <- c(.gettext("None (whole corpus)"), colnames(meta(corpus)))
    varBox <- variableListBox(top, vars,
                              title=.gettext("Report results by variable:"),
                              initialSelection=0)

    onOK <- function() {
        n <- as.numeric(tclvalue(tclN))
        termsList <- strsplit(tclvalue(tclTerms), " ")[[1]]
        var <- getSelection(varBox)

        if(length(termsList) == 0) {
            .Message(gettext("Please enter at least one term."), "error", parent=top)

            return()
        }
        else if(!all(termsList %in% colnames(dtm))) {
            wrongTerms <- termsList[!(termsList %in% colnames(dtm))]
            .Message(sprintf(.ngettext(length(wrongTerms),
                                      "Term \'%s\' does not exist in the corpus.",
                                      "Terms \'%s\' do not exist in the corpus."),
                                       # TRANSLATORS: this should be opening quote, comma, closing quote
                                       paste(wrongTerms, collapse=.gettext("\', \'"))),
                     "error", parent=top)

            return()
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        if(var == .gettext("None (whole corpus)")) {
            if(length(termsList) == 1)
                doItAndPrint(sprintf('coocs <- termChisqDist("%s", dtm, %i)', termsList, n))
            else
                doItAndPrint(sprintf('coocs <- sapply(c("%s"), termChisqDist, dtm, %i, simplify=FALSE)',
                                     paste(termsList, collapse='", "'), n))
        }
        else {
            if(length(termsList) == 1)
                doItAndPrint(sprintf('coocs <- termChisqDist("%s", dtm, %i, meta(corpus, "%s")[[1]])',
                                     termsList, n, var))
            else
                doItAndPrint(sprintf('coocs <- sapply(c("%s"), termChisqDist, dtm, %i, meta(corpus, "%s")[[1]], simplify=FALSE)',
                                     paste(termsList, collapse='", "'), n, var))
        }

        title <- sprintf(.ngettext(length(termsList),
                                   "Terms associated with term \"%s\" according to Chi-squared distance",
                                   "Terms associated with terms \"%s\" according to Chi-squared distance"),
                         # TRANSLATORS: this should be opening quote, comma, closing quote
                         paste(termsList, collapse=.gettext("\", \"")), n)

       if(var != .gettext("None (whole corpus)"))
           setLastTable("coocs", paste(title, sprintf(.gettext("(for %s)"),
                                     paste(levels(factor(meta(corpus, var)[[1]])),
                                           collapse=", "))))
       else
           setLastTable("coocs", title)

        doItAndPrint("coocs")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="termsCoocDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Reference terms (space-separated):")), sticky="w")
    tkgrid(entryTerms, sticky="w", columnspan=2)
    tkgrid(labelRcmdr(top, text=.gettext("Number of terms to show:")), spinN, sticky="sw", pady=6)
    tkgrid(getFrame(varBox), columnspan=2, sticky="w", pady=6)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix(focus=entryTerms)
}
