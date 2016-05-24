frequentTerms <- function(dtm, variable=NULL, n=25) {
    if(length(variable) == 1 && is.na(variable)) {
        counts <- sort(col_sums(dtm), decreasing=TRUE)[1:n]
        mat <- cbind(counts, counts/sum(dtm) * 100)
        colnames(mat) <- c(.gettext("Global freq."), .gettext("Global %"))
        return(mat)
    }

    if(!is.null(variable) && length(unique(variable)) < 2)
        stop("Please provide a variable with at least two levels.")

    if(!is.null(variable))
        dtm <- rollup(dtm, 1, variable)

    rs <- row_sums(dtm)
    tot <- sum(rs)

    cs <- col_sums(dtm)
    cs.tot <- cs/tot

    sapply(rownames(dtm), simplify=FALSE, function(l) {
        # rownames(dtm) == l is used below because "" is a possible level
        i <- rownames(dtm) == l

        # Empty documents create errors
        if(rs[i] == 0)
            return(numeric(0))

        rp <- as.matrix(dtm[l,]/rs[l])[1,]
        cp <- as.matrix(dtm[l,])[1,]/cs
        sup <- rp > cs.tot

        counts <- as.matrix(dtm[i,])[1,]

        # As this is a discrete distribution, we need to subtract one
        # to include the value when switching sides
        p.val <- phyper(ifelse(sup, counts - 1, counts), rs[l], tot - rs[l], cs)
        t.val <- qnorm(p.val)

        p.val[sup] <- 1 - p.val[sup]

        ord <- head(order(counts, decreasing=TRUE), n)
        ret <- cbind(term.clus=rp[ord] * 100, clus.term=cp[ord] * 100,
                     p.global=cs[ord]/tot * 100, n.int=counts[ord], n.global=cs[ord],
                     t.value=t.val[ord], p.value=round(p.val[ord], 4))
        colnames(ret) <- c(.gettext("% Term/Level"), .gettext("% Level/Term"), .gettext("Global %"),
                           .gettext("Level"), .gettext("Global"),
                           .gettext("t value"), .gettext("Prob."))
        ret[order(-sign(ret[, 6]), ret[, 7]),]
    })
}


freqTermsDlg <- function() {
    if(!(exists("dtm") && class(dtm) == "DocumentTermMatrix")) {
        .Message(message=.gettext("Please import a corpus and create the document-term matrix first."),
                type="error")
        return()
    }

    initializeDialog(title=.gettext("Show Most Frequent Terms"))
    tclN <- tclVar(10)
    spinN <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                      inc=1, textvariable=tclN,
                      validate="all", validatecommand=.validate.uint)

    vars <- c(.gettext("None (whole corpus)"), .gettext("Document"), colnames(meta(corpus)))
    varBox <- variableListBox(top, vars,
                              title=.gettext("Report results by variable:"),
                              initialSelection=0)

    onOK <- function() {
        var <- getSelection(varBox)
        n <- as.numeric(tclvalue(tclN))

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        if(var == .gettext("None (whole corpus)"))
            doItAndPrint(sprintf('freqTerms <- frequentTerms(dtm, NA, %i)', n))
        else if(var == .gettext("Document"))
            doItAndPrint(sprintf('freqTerms <- frequentTerms(dtm, NULL, %i)', n))
        else
            doItAndPrint(sprintf('freqTerms <- frequentTerms(dtm, meta(corpus, "%s")[[1]], %i)',
                                 var, n))

        if(var == .gettext("None (whole corpus)"))
            setLastTable("freqTerms", .gettext("Most frequent terms in the corpus"))
        else
            setLastTable("freqTerms", sprintf(.gettext("Most frequent terms by %s"), var))

        doItAndPrint("freqTerms")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="freqTermsDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Number of terms to show:")), spinN,
           sticky="sw", pady=6)
    tkgrid(getFrame(varBox), columnspan=2, sticky="w", pady=6)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix()
}
