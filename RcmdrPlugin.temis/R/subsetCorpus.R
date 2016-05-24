.subsetCorpus <- function(save) {
    if(save)
        doItAndPrint("origCorpus <- corpus")

    doItAndPrint("corpus <- corpus[keep]")

    if(exists("dtm")) {
        processing <- meta(corpus, type="corpus", tag="processing")
        lang <- meta(corpus, type="corpus", tag="language")

        if(save)
            doItAndPrint("origDtm <- dtm")

        doItAndPrint('dtmAttr <- attributes(dtm)')
        doItAndPrint('origDictionary <- attr(dtm, "dictionary")')

        doItAndPrint("dtm <- dtm[keep,]")

        # stemming=FALSE since we don't need to extract the stems again,
        # we reuse below those of the old dictionary
        .buildDictionary(FALSE, processing["customStemming"], lang)

        if(processing["stemming"] || processing["customStemming"])
            doItAndPrint('dictionary[[2]] <- origDictionary[rownames(dictionary), 2]')

        # customStemming=FALSE since we don't want to ask the user to customize stemming
        .prepareDtm(processing["stopwords"], processing["stemming"] || processing["customStemming"], FALSE, lang)
        if(processing["customStemming"])
            doItAndPrint('dtm <- dtm[, Terms(dtm) != ""]')

        doItAndPrint('attr(dtm, "language") <- dtmAttr$lang')
        doItAndPrint('attr(dtm, "processing") <- dtmAttr$processing')

        doItAndPrint("rm(dtmAttr, origDictionary)")
    }

    if(exists("wordsDtm")) {
        if(save)
            doItAndPrint("origWordsDtm <- wordsDtm")

        doItAndPrint("wordsDtm <- wordsDtm[keep,]")
        doItAndPrint("wordsDtm <- wordsDtm[,col_sums(wordsDtm) > 0]")
    }

    doItAndPrint("corpusVars <- corpusVars[keep,, drop=FALSE]")

    # Remove objects left from a previous analysis on the old corpus to avoid confusion
    # (we assume later existing objects match the current corpus)
    objects <- c("keep", "voc", "lengths", "termFreqs", "corpusClust", "corpusSubClust", "corpusCa", "plottingCa")
    doItAndPrint(paste('rm(list=c("', paste(objects[sapply(objects, exists)], collapse='", "'), '"))', sep=""))
    gc()

    doItAndPrint("corpus")
    doItAndPrint("dtm")
}

subsetCorpusByVarDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Subset Corpus by Variable"))

    vars <- colnames(meta(corpus))

    # We cannot use variableListBox as it is not meant for changing levels
    varsFrame <- tkframe(top)
    varsBox <- tklistbox(varsFrame, height=getRcmdr("variable.list.height"),
                         selectmode="single", export=FALSE)
    varsScrollbar <- ttkscrollbar(varsFrame, command=function(...) tkyview(varsBox, ...))
    tkconfigure(varsBox, yscrollcommand=function(...) tkset(varsScrollbar, ...))
    for(var in vars) tkinsert(varsBox, "end", var)
    tkselection.set(varsBox, 0)

    levelsFrame <- tkframe(top)
    levelsBox <- tklistbox(levelsFrame, height=getRcmdr("variable.list.height"),
                           selectmode=getRcmdr("multiple.select.mode"), export=FALSE)
    levelsScrollbar <- ttkscrollbar(levelsFrame, command=function(...) tkyview(levelsBox, ...))
    tkconfigure(levelsBox, yscrollcommand=function(...) tkset(levelsScrollbar, ...))
    for(level in unique(meta(corpus, vars[1])[[1]])) tkinsert(levelsBox, "end", level)
    tkselection.set(levelsBox, 0)

    onSelect <- function() {
        var <- vars[as.numeric(tkcurselection(varsBox))+1]
        tkdelete(levelsBox, "0", "end")

        levs <- unique(meta(corpus, var)[[1]])
        for(level in levs) tkinsert(levelsBox, "end", level)
    }

    tkbind(varsBox, "<<ListboxSelect>>", onSelect)

    tclSave <- tclVar("1")
    checkSave <- tkcheckbutton(top, text=.gettext("Save original corpus to restore it later"),
                               variable=tclSave)

    onOK <- function() {
        var <- vars[as.numeric(tkcurselection(varsBox))+1]
        levs <- unique(meta(corpus, var)[[1]])[as.numeric(tkcurselection(levelsBox))+1]
        save <- tclvalue(tclSave) == "1"

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        doItAndPrint(sprintf('keep <- meta(corpus, "%s")[[1]] %%in%% c("%s")',
                             var, paste(levs, collapse='", "')))

        .subsetCorpus(save)

	    activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="subsetCorpusByVarDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Select a variable and one or more levels to retain:")),
           columnspan=2, sticky="w", pady=6)
    tkgrid(.titleLabel(varsFrame, text=.gettext("Variable:")),
           sticky="w")
    tkgrid(varsBox, varsScrollbar, sticky="ewns", pady=6)
    tkgrid(.titleLabel(levelsFrame, text=.gettext("Levels:")),
           sticky="w")
    tkgrid(levelsBox, levelsScrollbar, sticky="ewns", pady=6)
    tkgrid(varsFrame, levelsFrame, sticky="wns", pady=6)
    tkgrid(checkSave, sticky="w", pady=6)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix(focus=varsBox)
}


subsetCorpusByTermsDlg <- function() {
    initializeDialog(title=.gettext("Subset Corpus by Terms"))

    tclKeep <- tclVar("")
    entryKeep <- ttkentry(top, width="40", textvariable=tclKeep)

    tclKeepFreq <- tclVar(1)
    spinKeep <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                         inc=1, textvariable=tclKeepFreq,
                         validate="all", validatecommand=.validate.uint)

    tclExclude <- tclVar("")
    entryExclude <- ttkentry(top, width="40", textvariable=tclExclude)

    tclExcludeFreq <- tclVar(1)
    spinExclude <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                            inc=1, textvariable=tclExcludeFreq,
                            validate="all", validatecommand=.validate.uint)

    tclSave <- tclVar("1")
    checkSave <- tkcheckbutton(top, text=.gettext("Save original corpus to restore it later"),
                               variable=tclSave)

    onOK <- function() {
        keepList <- strsplit(tclvalue(tclKeep), " ")[[1]]
        excludeList <- strsplit(tclvalue(tclExclude), " ")[[1]]
        keepFreq <- as.numeric(tclvalue(tclKeepFreq))
        excludeFreq <- as.numeric(tclvalue(tclExcludeFreq))
        save <- tclvalue(tclSave) == "1"

        if(length(keepList) == 0 && length(excludeList) == 0) {
            .Message(.gettext("Please enter at least one term."), "error", parent=top)

            return()
        }
        else if(!all(c(keepList, excludeList) %in% colnames(dtm))) {
            wrongTerms <- c(keepList, excludeList)[!c(keepList, excludeList) %in% colnames(dtm)]
            .Message(sprintf(.ngettext(length(wrongTerms),
                                      "Term \'%s\' does not exist in the corpus.",
                                      "Terms \'%s\' do not exist in the corpus."),
                                      # TRANSLATORS: this should be opening quote, comma, closing quote
                                      paste(wrongTerms, collapse=.gettext("\', \'"))),
                     "error", parent=top)

            return()
        }
        else if((length(keepList) > 0 && length(excludeList) > 0 &&
                 !any(row_sums(dtm[, keepList] >= keepFreq) > 0 &
                      row_sums(dtm[, excludeList] >= excludeFreq) == 0)) ||
                (length(keepList) > 0 && length(excludeList) == 0 &&
                 !any(row_sums(dtm[, keepList] >= keepFreq) > 0)) ||
                (length(keepList) == 0 && length(excludeList) > 0 &&
                 !any(row_sums(dtm[, excludeList] >= excludeFreq) == 0))) {
            .Message(.gettext("Specified conditions would exclude all documents from the corpus."),
                     "error", parent=top)

            return()
        }

        closeDialog()

        if(length(keepList) > 0 && length(excludeList) > 0)
            doItAndPrint(sprintf('keep <- row_sums(dtm[, c("%s")] >= %i) > 0 & row_sums(dtm[, c("%s")] >= %i) == 0',
                                 paste(keepList, collapse='", "'), keepFreq,
                                 paste(excludeList, collapse='", "'), excludeFreq))
        else if(length(keepList) > 0)
            doItAndPrint(sprintf('keep <- row_sums(dtm[, c("%s")] >= %i) > 0',
                                 paste(keepList, collapse='", "'), keepFreq))
        else
            doItAndPrint(sprintf('keep <- row_sums(dtm[, c("%s")] >= %i) == 0',
                                 paste(excludeList, collapse='", "'), excludeFreq))

        .subsetCorpus(save)

	    activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="subsetCorpusByTermsDlg")
    tkgrid(labelRcmdr(top, text=.gettext("Keep documents containing one of these terms (space-separated):")),
           sticky="w", columnspan=4)
    tkgrid(entryKeep,
           labelRcmdr(top, text=.gettext("at least")),
           spinKeep,
           labelRcmdr(top, text=.gettext("time(s)")),
           sticky="w", pady=c(0, 6))
    tkgrid(labelRcmdr(top, text=.gettext("Exclude documents containing one of these terms (space-separated):")),
           sticky="w", pady=c(6, 0), columnspan=4)
    tkgrid(entryExclude,
           labelRcmdr(top, text=.gettext("at least")),
           spinExclude,
           labelRcmdr(top, text=.gettext("time(s)")),
           sticky="w", pady=c(0, 6))
    tkgrid(labelRcmdr(top, text=.gettext("(Only documents matching both conditions will be retained in the new corpus.)")),
           sticky="w", pady=6, columnspan=4)
    tkgrid(checkSave, sticky="w", pady=c(12, 6), columnspan=4)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=4)
    dialogSuffix(focus=entryKeep)
}

restoreCorpus <- function() {
    if(!exists("origCorpus"))
        .Message(message=.gettext("No saved corpus to restore was found."), type="error")

    doItAndPrint("corpus <- origCorpus")

    if(exists("origDtm"))
        doItAndPrint("dtm <- origDtm")

    if(exists("origWordsDtm"))
        doItAndPrint("dtm <- origWordsDtm")


    # Remove objects left from a previous analysis on the subset corpus to avoid confusion
    # (we assume later existing objects match the current corpus)
    objects <- c("keep", "voc", "lengths", "termFreqs", "absTermFreqs", "varTermFreqs",
                 "corpusClust", "corpusSubClust", "corpusCa", "plottingCa",
                 # Also remove backup objects
                 "origCorpus", "origDtm", "origWordsDtm")
    doItAndPrint(paste("rm(", paste(objects[sapply(objects, exists)], collapse=", "), ")", sep=""))

    doItAndPrint("corpus")
    doItAndPrint("dtm")
    gc()
}

