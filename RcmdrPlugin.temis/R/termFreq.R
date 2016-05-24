termFrequencies <- function(dtm, terms, variable=NULL, n=25, by.term=FALSE) {
   wrongTerms <- terms[!terms %in% colnames(dtm)]
   if(length(wrongTerms) > 0)
       stop(sprintf(.ngettext(length(wrongTerms),
                              "Term \'%s\' does not exist in the corpus.",
                              "Terms \'%s\' do not exist in the corpus."),
                    paste(wrongTerms, collapse=.gettext("\', \'"))))

    if(length(variable) == 1 && is.na(variable)) {
        counts <- col_sums(dtm[, terms])
        mat <- cbind(counts, counts/sum(dtm) * 100)
        colnames(mat) <- c(.gettext("Global"), .gettext("Global %"))
        return(mat)
    }

    if(!is.null(variable) && length(unique(variable)) < 2)
        stop("Please provide a variable with at least two levels.")

    if(!is.null(variable))
        dtm <- rollup(dtm, 1, variable)

    # We need to compute these statistics before removing terms so that they are stable
    rs <- row_sums(dtm)
    tot <- sum(rs)

    dtm <- dtm[, terms]

    cs <- col_sums(dtm)
    cs.tot <- cs/tot

    ret <- sapply(rownames(dtm), simplify="array", function(l) {
        # rownames(dtm) == l is used below because "" is a possible level
        i <- rownames(dtm) == l

        rp <- as.matrix(dtm[l,]/rs[l])[1,]
        cp <- as.matrix(dtm[l,])[1,]/cs
        sup <- rp > cs.tot

        counts <- as.matrix(dtm[i,])[1,]

        # As this is a discrete distribution, we need to subtract one
        # to include the value when switching sides
        p.val <- phyper(ifelse(sup, counts - 1, counts), rs[l], tot - rs[l], cs)
        t.val <- qnorm(p.val)

        p.val[sup] <- 1 - p.val[sup]

        ret <- cbind(term.clus=rp * 100, clus.term=cp * 100,
                     p.global=cs/tot * 100, n.int=counts, n.global=cs,
                     t.value=t.val, p.value=round(p.val, 4))
        colnames(ret) <- c(.gettext("% Term/Level"), .gettext("% Level/Term"), .gettext("Global %"),
                           .gettext("Level"), .gettext("Global"),
                           .gettext("t value"), .gettext("Prob."))
        ret
    })

   if(length(terms) == 1 || by.term == TRUE)
       aperm(ret, c(3, 2, 1))
   else
       ret
}

termFreqDlg <- function() {
    initializeDialog(title=.gettext("Term Frequencies"))

    tclTerms <- tclVar("")
    entryTerms <- ttkentry(top, width="30", textvariable=tclTerms)

    vars <- c(.gettext("Document"), colnames(meta(corpus)))
    varBox <- variableListBox(top, vars,
                              title=.gettext("Variable:"),
                              initialSelection=0)

    radioButtons(name="what",
                 buttons=c("term.lev", "lev.term", "occ"),
                 labels=c(.gettext("Term prevalence in level (\"% Term/Level\")"),
                          .gettext("Distribution of occurrences among levels (\"% Level/Term\")"),
                          .gettext("Absolute number of occurrences in level (\"Level\")")),
                 title=.gettext("Measure to plot:"),
                 right.buttons=FALSE)

    displayFrame <- tkframe(top)

    tclTitle <- tclVar(.gettext("Occurrences of term %T by %V"))
    titleEntry <- ttkentry(displayFrame, width="40", textvariable=tclTitle)


    tclByTermVar <- tclVar(1)
    transButton <- tkcheckbutton(displayFrame, text=.gettext("Group results by term"), variable=tclByTermVar)

    tclPlotVar <- tclVar(1)
    plotButton <- tkcheckbutton(displayFrame, text=.gettext("Draw plot"), variable=tclPlotVar)

    tclVertVar <- tclVar(0)
    vertButton <- tkcheckbutton(displayFrame, text=.gettext("Vertical bars"), variable=tclVertVar)

    tclStackVar <- tclVar(0)
    stackButton <- tkcheckbutton(displayFrame, text=.gettext("Stacked bars"), variable=tclStackVar)

    onOK <- function() {
        termsList <- strsplit(tclvalue(tclTerms), " ")[[1]]
        var <- getSelection(varBox)
        title <- tclvalue(tclTitle)
        plot <- tclvalue(tclPlotVar) == 1
        vert <- tclvalue(tclVertVar) == 1
        stack <- tclvalue(tclStackVar) == 1
        byTerm <- tclvalue(tclByTermVar) == 1

        if(length(termsList) == 0) {
            .Message(.gettext("Please enter at least one term."), "error", parent=top)

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

        what <- tclvalue(whatVariable)

        # Table
        if(!byTerm) {
            if(var == .gettext("None (whole corpus)"))
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), NA)',
                                     paste(termsList, collapse='", "')))
            else if(var == .gettext("Document"))
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), NULL)',
                                     paste(termsList, collapse='", "')))
            else
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), meta(corpus, "%s")[[1]])',
                                     paste(termsList, collapse='", "'), var))
       }
       else {
            if(var == .gettext("None (whole corpus)"))
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), NA, by.term=TRUE)',
                                     paste(termsList, collapse='", "')))
            else if(var == .gettext("Document"))
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), NULL, by.term=TRUE)',
                                     paste(termsList, collapse='", "')))
            else
                doItAndPrint(sprintf('termFreqs <- termFrequencies(dtm, c("%s"), meta(corpus, "%s")[[1]], by.term=TRUE)',
                                     paste(termsList, collapse='", "'), var))
       }


        if(what == "term.lev")
            lab <- .gettext("% of all terms")
        else if (what == "lev.term")
            lab <- .gettext("% of occurrences")
        else
            lab <- .gettext("Number of occurrences")

        if(length(termsList) == 1)
            title <- gsub("%T", termsList[1], title)
        else
            title <- gsub(" %T ", " ", title)

        title <- gsub("%V", tolower(var), title)

        # Plot
        if(plot) {
            col <- switch(what, term.lev=.gettext("% Term/Level"), lev.term=.gettext("% Level/Term"), .gettext("Level"))

            if(vert) {
                if(length(termsList) > 1)
                    doItAndPrint(sprintf('barchart(t(termFreqs[, "%s",]), horizontal=%s, stack=%s, %s="%s", main="%s", auto.key=TRUE, scales=list(rot=90))',
                                         col, !vert, stack, if(vert) "ylab" else "xlab", lab, title))
                else
                    doItAndPrint(sprintf('barchart(cbind(termFreqs[, "%s",]), horizontal=%s, %s="%s", main="%s", scales=list(rot=90))',
                                         col, !vert, if(vert) "ylab" else "xlab", lab, title))
            }
            else {
                if(length(termsList) > 1)
                    doItAndPrint(sprintf('barchart(t(termFreqs[, "%s",]), stack=%s, xlab="%s", main="%s", auto.key=TRUE)',
                                         col, stack, lab, title))
                else
                    doItAndPrint(sprintf('barchart(termFreqs[, "%s",], xlab="%s", main="%s")',
                                         col, lab, title))
            }
        }

        if(what == "term.lev")
            setLastTable("termFreqs", paste(title, .gettext("(% of all terms)")))
        else if(what == "lev.term")
            setLastTable("termFreqs", paste(title, .gettext("(% of occurrences)")))
        else
            setLastTable("termFreqs", title)

        doItAndPrint("termFreqs")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="termFreqDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Terms to show (space-separated):")), sticky="w", columnspan=2)
    tkgrid(entryTerms, sticky="w", columnspan=2)
    tkgrid(getFrame(varBox), sticky="w", columnspan=2, pady=6)
    tkgrid(whatFrame, sticky="w", columnspan=2, pady=6)
    tkgrid(.titleLabel(displayFrame, text=.gettext("Display:")),
           sticky="w", columnspan=2)
    tkgrid(labelRcmdr(displayFrame, text=.gettext("Title:")), titleEntry, sticky="w", padx=6)
    tkgrid(transButton, sticky="w", columnspan=2)
    tkgrid(plotButton, sticky="w", columnspan=2)
    tkgrid(vertButton, sticky="w", columnspan=2)
    tkgrid(stackButton, sticky="w", columnspan=2)
    tkgrid(displayFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=entryTerms)
}
