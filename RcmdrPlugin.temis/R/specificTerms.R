specificTerms <- function(dtm, variable=NULL, p=0.1, n.max=25, sparsity=0.95, min.occ=2) {
    if(!is.null(variable) && length(unique(variable)) < 2)
        stop(.gettext("Please provide a variable with at least two levels."))

    if(!is.null(variable))
        dtm <- rollup(dtm, 1, variable)

    # We need to compute these statistics before removing terms so that they are stable
    rs <- row_sums(dtm)
    tot <- sum(rs)

    if(sparsity < 1)
        dtm <- removeSparseTerms(dtm, sparsity)

    if(min.occ > 1)
        dtm <- dtm[, col_sums(dtm) >= min.occ]

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

        keep <- which(p.val <= p)

        if(length(keep) == 0) return(numeric(0))

        ord <- head(intersect(order(p.val), keep), n.max)
        ret <- cbind(term.clus=rp[ord] * 100, clus.term=cp[ord] * 100,
                     p.global=cs[ord]/tot * 100, n.int=counts[ord], n.global=cs[ord],
                     t.value=t.val[ord], p.value=round(p.val[ord], 4))
        colnames(ret) <- c(.gettext("% Term/Level"), .gettext("% Level/Term"), .gettext("Global %"),
                           .gettext("Level"), .gettext("Global"),
                           .gettext("t value"), .gettext("Prob."))

        ret[order(-sign(ret[, 6]), ret[, 7]), , drop=FALSE]
    })
}

specificTermsDlg <- function() {
    if(!(exists("dtm") && class(dtm) == "DocumentTermMatrix")) {
        .Message(message=.gettext("Please import a corpus and create the document-term matrix first."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Show Specific Terms"))

    vars <- c(.gettext("Document"), colnames(meta(corpus)))
    varBox <- variableListBox(top, vars,
                              title=.gettext("Show terms specific of levels of variable:"),
                              initialSelection=0)

    tclP <- tclVar(10)
    spinP <- tkwidget(top, type="spinbox", from=0, to=100,
                      inc=1, textvariable=tclP,
                      validate="all", validatecommand=.validate.unum)

    tclOcc <- tclVar(2)
    spinOcc <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                        inc=1, textvariable=tclOcc,
                        validate="all", validatecommand=.validate.uint)

    tclN <- tclVar(25)
    spinN <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                      inc=1, textvariable=tclN,
                      validate="all", validatecommand=.validate.uint)

    onOK <- function() {
        var <- getSelection(varBox)
        p <- as.numeric(tclvalue(tclP))
        occ <- as.numeric(tclvalue(tclOcc))
        n <- as.numeric(tclvalue(tclN))

        if(var != .gettext("Document") && length(unique(meta(corpus, var)[[1]])) < 2) {
            .Message(.gettext("Please choose a variable with at least two levels."), "error", parent=top)
            return()
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        if(var == .gettext("Document")) {
            doItAndPrint(sprintf('specTerms <- specificTerms(dtm, NULL, p=%s, min.occ=%s, n.max=%s)',
                                 p/100, occ, n))
        }
        else {
            doItAndPrint(sprintf('specTerms <- specificTerms(dtm, meta(corpus, "%s")[[1]], p=%s, min.occ=%s, n.max=%s)',
                                 var, p/100, occ, n))
        }

        if(var == .gettext("Document"))
            setLastTable("specTerms", .gettext("Specific terms by document"))
        else
            setLastTable("specTerms", sprintf(.gettext("Specific terms by %s"), var))

        doItAndPrint("specTerms")

        activateMenus()

        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="specificTermsDlg")
    tkgrid(getFrame(varBox), columnspan=2, sticky="w", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Show terms with a probability below (%):")), spinP,
           sticky="sw", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Only retain terms with a number of occurrences above:")), spinOcc,
           sticky="sw", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Maximum number of terms to show per level:")), spinN,
           sticky="sw", pady=6)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix()
}

