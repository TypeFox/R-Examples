recodeTimeVarDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Recode Time Variable"))

    vars <- colnames(meta(corpus))[colnames(meta(corpus)) != "MetaID"]

    datevar <- which(vars == .gettext("Date")) - 1
    timevar <- which(vars == .gettext("Time")) - 1
    datetimevar <- if(length(datevar) > 0) datevar
                   else if(length(timevar) > 0) timevar
                   else 0

    timeVarBox <- variableListBox(top, vars,
                                  title=.gettext("Existing date/time variable:"),
                                  initialSelection=datetimevar)

    inFrame <- tkframe(top)
    tclInFormat <- if(length(timevar) == 0) tclVar("%Y-%m-%d") else tclVar("%Y-%m-%d %H:%M")
    inFormatEntry <- ttkentry(inFrame, width="20", textvariable=tclInFormat)

    tclOutFormat <- if(length(timevar) == 0) tclVar("%Y-%m-%d") else tclVar("%Y-%m-%d %H:%M")
    outFormatEntry <- ttkentry(top, width="20", textvariable=tclOutFormat)

    tclNewName <- tclVar(paste(getSelection(timeVarBox), "2", sep=""))
    newNameEntry <- ttkentry(top, textvariable=tclNewName)


    onSelectTimeVar <- function() {
        var <- getSelection(timeVarBox)

        if(var == .gettext("Date"))
            format <- "%Y-%m-%d"
        else
            format <- "%Y-%m-%d %H:%M"

        tkdelete(inFormatEntry, "0", "end")
        tkinsert(inFormatEntry, "end", format)

        tkdelete(outFormatEntry, "0", "end")
        tkinsert(outFormatEntry, "end", format)

        tkdelete(newNameEntry, "0", "end")
        tkinsert(newNameEntry, "end", paste(var, "2", sep=""))
    }

    tkbind(timeVarBox$listbox, "<<ListboxSelect>>", onSelectTimeVar)

    onOK <- function() {
        timeVar <- getSelection(timeVarBox)
        newName <- tclvalue(tclNewName)
        inFormat <- tclvalue(tclInFormat)
        outFormat <- tclvalue(tclOutFormat)

        # Check that input format is more or less correct before running the code
        time <- meta(corpus, timeVar)[[1]]
        time <- strptime(unique(time[!is.na(time)]), inFormat)
        if(all(is.na(time))) {
            .Message(message=sprintf(.gettext("Incorrect input time format or variable: no values of \"%s\" could be converted to a time index."), timeVar),
                     type="error", parent=top)
            return()
        }
        else if(any(is.na(time))) {
            .Message(message=sprintf(.gettext("Some values of \"%s\" could not be converted to a time index and will be missing."), timeVar),
                     type="warning", parent=top)
        }

        closeDialog()

        doItAndPrint(sprintf('meta(corpus, "%s") <- strftime(strptime(meta(corpus, "%s")[[1]], format="%s"), format="%s")',
                             newName, timeVar, inFormat, outFormat))
        doItAndPrint(sprintf('corpusVars[["%s"]] <- meta(corpus, "%s")', newName, newName))
        doItAndPrint(sprintf('str(meta(corpus, "%s")[[1]])', newName))

        activateMenus()
        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="recodeTimeVarDlg")
    tkgrid(getFrame(timeVarBox), sticky="ewns", pady=6, padx=c(0, 6), row=0, rowspan=2)
    tkgrid(labelRcmdr(inFrame, text=.gettext("Time format of existing variable:")),
           pady=c(6, 0), sticky="ws")
    tkgrid(inFormatEntry, sticky="ws")
    tkgrid(inFrame, sticky="wn", padx=c(6, 0), row=0, column=1, rowspan=2)
    tkgrid(.titleLabel(top, text=.gettext("New variable name:")),
           sticky="ewns", padx=c(0, 6), pady=c(6, 0), row=2, column=0)
    tkgrid(newNameEntry, sticky="ewns", padx=c(0, 6), pady=c(0, 6), row=3, column=0)
    tkgrid(labelRcmdr(top, text=.gettext("Time format of new variable:")),
           padx=c(6, 0), pady=c(6, 0), sticky="w", row=2, column=1)
    tkgrid(outFormatEntry, sticky="w", padx=c(6, 0), pady=c(0, 6), row=3, column=1)
    tkgrid(labelRcmdr(top, text=.gettext("Useful codes:\n%Y: year - %m: month number - %B: month name\n%W: week number starting on Mondays - %U: starting on Sundays\n%d: day number - %A: week day name\n%H: hour - %M: minute\n\nClick the  \"Help\" button for more codes.")),
           pady=6, sticky="ewns", columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=timeVarBox$listbox)
}


varTimeSeriesDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Corpus Temporal Evolution"))

    vars <- colnames(meta(corpus))[colnames(meta(corpus)) != "MetaID"]

    datevar <- which(vars == .gettext("Date")) - 1
    timevar <- which(vars == .gettext("Time")) - 1
    datetimevar <- if(length(datevar) > 0) datevar
                   else if(length(timevar) > 0) timevar
                   else 0

    timeVarBox <- variableListBox(top, vars,
                                  title=.gettext("Date/time variable:"),
                                  initialSelection=datetimevar)

    tclFormat <- if(length(timevar) == 0) tclVar("%Y-%m-%d") else tclVar("%Y-%m-%d %H:%M")
    formatEntry <- ttkentry(top, width="20", textvariable=tclFormat)

    # We cannot use variableListBox as it is not meant for changing levels
    varsFrame <- tkframe(top)
    varsBox <- tklistbox(varsFrame, height=getRcmdr("variable.list.height"),
                         selectmode="single", export=FALSE)
    varsScrollbar <- ttkscrollbar(varsFrame, command=function(...) tkyview(varsBox, ...))
    tkconfigure(varsBox, yscrollcommand=function(...) tkset(varsScrollbar, ...))
    for(var in c(.gettext("None (one curve)"), vars)) tkinsert(varsBox, "end", var)
    tkselection.set(varsBox, 0)

    levelsFrame <- tkframe(top)
    levelsBox <- tklistbox(levelsFrame, height=getRcmdr("variable.list.height"),
                           selectmode=getRcmdr("multiple.select.mode"), export=FALSE)
    levelsScrollbar <- ttkscrollbar(levelsFrame, command=function(...) tkyview(levelsBox, ...))
    tkconfigure(levelsBox, yscrollcommand=function(...) tkset(levelsScrollbar, ...))
    tkselection.set(levelsBox, 0)

    onSelectTimeVar <- function() {
        var <- getSelection(timeVarBox)

        if(var == .gettext("Date")) {
            tkdelete(formatEntry, "0", "end")
            tkinsert(formatEntry, "end", "%Y-%m-%d")
        }
        else if (var == .gettext("Time")) {
            tkdelete(formatEntry, "0", "end")
            tkinsert(formatEntry, "end", "%Y-%m-%d %H:%M")
        }
    }

    tkbind(timeVarBox$listbox, "<<ListboxSelect>>", onSelectTimeVar)

    onSelectGroup <- function() {
        var <- c("", vars)[as.numeric(tkcurselection(varsBox))+1]
        tkdelete(levelsBox, "0", "end")

        if(var == "")
            return()

        levs <- if(is.factor(meta(corpus, var)[[1]])) levels(droplevels(meta(corpus, var)[[1]]))
                else sort(unique(meta(corpus, var)[[1]]))
        for(level in as.character(levs[!is.na(levs)])) tkinsert(levelsBox, "end", level)

        tkselection.set(levelsBox, 0, "end")
    }

    tkbind(varsBox, "<<ListboxSelect>>", onSelectGroup)

    radioButtons(name="what",
                 buttons=c("number", "percent"),
                 labels=c(.gettext("Number of documents per time unit"),
                          .gettext("% of documents per time unit")),
                 title=.gettext("Measure:"),
                 right.buttons=FALSE)

    tclMean <- tclVar(0)
    meanButton <- tkcheckbutton(top, text=.gettext("Apply rolling mean"), variable=tclMean)

    tclWindow <- tclVar(7)
    spinWindow <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                             inc=1, textvariable=tclWindow,
                             validate="all", validatecommand=.validate.uint)

    tclTitle <- tclVar(.gettext("Temporal evolution of the corpus"))
    titleEntry <- ttkentry(top, width="30", textvariable=tclTitle)

    onCustom <- function() {
        timeVar <- getSelection(timeVarBox)
        groupVar <- c("", vars)[as.numeric(tkcurselection(varsBox))+1]
        what <- tclvalue(whatVariable)
        rollmean <- tclvalue(tclMean) == 1
        window <- as.numeric(tclvalue(tclWindow))
        title <- tclvalue(tclTitle)

        if(what == "percent" && nchar(groupVar) == 0) {
            .Message(message=.gettext("Plotting percents of documents with only one curve does not make sense: all points would be 100%."),
                     type="error", parent=top)
            return()
        }

        if(nchar(groupVar) > 0)
            groupLevs <- unique(meta(corpus, groupVar)[[1]])[as.numeric(tkcurselection(levelsBox))+1]

        format <- tclvalue(tclFormat)

        # Check that format is more or less correct before running the code
        time <- meta(corpus, timeVar)[[1]]
        time <- strptime(unique(time[!is.na(time)]), format)
        if(all(is.na(time))) {
            .Message(message=sprintf(.gettext("Incorrect time format or variable: no values of \"%s\" could be converted to a time index."), timeVar),
                     type="error", parent=top)
            return()
        }
        else if(any(is.na(time))) {
            .Message(message=sprintf(.gettext("Some values of \"%s\" could not be converted to a time index and will be missing."), timeVar),
                     type="warning", parent=top)
        }

        if(nchar(groupVar) == 0) {
            doItAndPrint(sprintf('tab <- table(as.character(strptime(meta(corpus, "%s")[[1]], "%s")))', timeVar, format))
            doItAndPrint("time <- as.POSIXct(names(tab))")
            doItAndPrint("docSeries <- zoo(tab, order.by=time)")
        }
        else {
            doItAndPrint(sprintf('tab <- table(as.character(strptime(meta(corpus, "%s")[[1]], "%s")), meta(corpus, "%s")[[1]])',
                                               timeVar, format, groupVar))

            if(what == "percent")
                doItAndPrint("tab <- prop.table(tab, 1)*100")

            doItAndPrint("time <- as.POSIXct(rownames(tab))")

            if(length(groupLevs) < length(unique(meta(corpus, var)[[1]])))
                doItAndPrint(sprintf('docSeries <- zoo(tab[,c("%s")], order.by=time)',
                                     paste(groupLevs, collapse='", "')))
            else
                doItAndPrint("docSeries <- zoo(tab, order.by=time)")
        }


        # We need to be sure we have a common time unit, i.e. we have a zooreg object
        if(length(docSeries) > 1 && !is.regular(docSeries))
            doItAndPrint('docSeries <- aggregate(docSeries, list(as.POSIXct(trunc(time, units(diff(time))))), regular=TRUE)')


        # For some reason, computing this after merging returns 24 "hours" instead of 1 "day" as unit
        unitsLoc <- c(secs=.gettext("per second"), mins=.gettext("per minute"), hours=.gettext("per hour"),
                      days=.gettext("per day"), weeks=.gettext("per week"))
        unit <- unitsLoc[units(diff(time(docSeries)))]

        # Trick to get a strictly regular time series with 0 where no document was observed
        # difftime chooses the unit so that all differences are > 1, which is what we want
        if(length(docSeries) > 1 && !is.regular(docSeries, strict=TRUE)) {
            # seq() will specify different if we pass "day"
            byUnit <- units(diff(time(docSeries)))
            if(byUnit == "days") byUnit <- "DSTday"

            doItAndPrint(sprintf('docSeries <- merge(docSeries, zoo(, seq(start(docSeries), end(docSeries), "%s")), fill=%s)',
                                 byUnit, if(what == "number") "0" else "NaN"))
        }

        if(rollmean) {
            if(window >= NROW(docSeries))
                .Message(message=.gettext("Chosen roll mean window is longer than the range of the time variable, rolling mean was not applied."),
                        type="warning", parent=top)
            else
                # For percents, the days with no observation get 0/0 == NaN, and we need to skip them
                doItAndPrint(sprintf('docSeries <- rollapply(docSeries, %s, align="left", mean, na.rm=TRUE)', window))
        }

        ylab <- if(what == "number") .gettext("Number of documents") else .gettext("% of documents")
        if(length(docSeries) > 1)
            doItAndPrint(sprintf('xyplot(docSeries, superpose=TRUE, xlab="", ylab="%s", main="%s", auto.key=%s, par.settings=simpleTheme(lwd=1.5))',
                                 paste(ylab, unit), title,
                                 if(NCOL(docSeries) > 1) 'TRUE' else "NULL"))
        else
            .Message(.gettext("Only one time point present, no plot can be drawn."), "error", parent=top)

        doItAndPrint("rm(tab, time)")

        setLastTable("docSeries", paste(title, " (", ylab, " ", unit, ")", sep=""))

        doItAndPrint("docSeries")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    # Shut up R CMD check WARNING
    onClose <- NULL
    .customCloseHelp(helpSubject="varTimeSeriesDlg", custom.button=.gettext("Draw plot"))

    tkgrid(getFrame(timeVarBox), sticky="ewns", padx=c(0, 6), pady=6, row=0, rowspan=3)
    tkgrid(labelRcmdr(top, text=.gettext("Time format:")), padx=c(6, 0), pady=c(6, 0), sticky="w", row=0, column=1)
    tkgrid(formatEntry, sticky="w", padx=c(6, 0), row=1, column=1)
    tkgrid(labelRcmdr(top, text=.gettext("%Y: year - %m: month - %d: day\n%H: hour - %M: minute\nClick the \"Help\" button for more codes.")),
           sticky="w", row=2, column=1, padx=c(6, 0), pady=6)
    tkgrid(.titleLabel(varsFrame, text=.gettext("Group by variable:")),
           sticky="w", pady=c(12, 0))
    tkgrid(varsBox, varsScrollbar, sticky="ewns", pady=6)
    tkgrid(varsFrame, levelsFrame, sticky="ewns", pady=6)
    tkgrid(labelRcmdr(levelsFrame, text=.gettext("Only plot levels:")), sticky="w", pady=c(12, 0))
    tkgrid(levelsBox, levelsScrollbar, sticky="ewns", pady=6)
    tkgrid(whatFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(.titleLabel(top, text=.gettext("Rolling mean:")),
           sticky="w", pady=c(6, 0))
    tkgrid(meanButton, sticky="w")
    tkgrid(labelRcmdr(top, text=.gettext("Time window for mean (in time units):")), spinWindow, sticky="w",
           padx=6, pady=c(0, 6))
    tkgrid(.titleLabel(top, text=.gettext("Title:")),
           sticky="w", pady=c(6, 0))
    tkgrid(titleEntry, sticky="w", padx=6, pady=c(0, 6), columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=timeVarBox$listbox, onOK=onCustom, onCancel=onClose)
}

termTimeSeriesDlg <- function() {
    nVars <- ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"])
    if(nVars == 0) {
        .Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Temporal Evolution of Term Occurrences"))

    vars <- colnames(meta(corpus))[colnames(meta(corpus)) != "MetaID"]

    datevar <- which(vars == .gettext("Date")) - 1
    timevar <- which(vars == .gettext("Time")) - 1
    datetimevar <- if(length(datevar) > 0) datevar
                   else if(length(timevar) > 0) timevar
                   else 0

    timeVarBox <- variableListBox(top, vars,
                                  title=.gettext("Date/time variable:"),
                                  initialSelection=datetimevar)

    tclFormat <- if(length(timevar) == 0) tclVar("%Y-%m-%d") else tclVar("%Y-%m-%d %H:%M")
    formatEntry <- ttkentry(top, width="20", textvariable=tclFormat)

    # We cannot use variableListBox as it is not meant for changing levels
    varsFrame <- tkframe(top)
    varsBox <- tklistbox(varsFrame, height=getRcmdr("variable.list.height"),
                         selectmode="single", export=FALSE)
    varsScrollbar <- ttkscrollbar(varsFrame, command=function(...) tkyview(varsBox, ...))
    tkconfigure(varsBox, yscrollcommand=function(...) tkset(varsScrollbar, ...))
    for(var in c(.gettext("None (one curve)"), vars)) tkinsert(varsBox, "end", var)
    tkselection.set(varsBox, 0)

    levelsFrame <- tkframe(top)
    levelsBox <- tklistbox(levelsFrame, height=getRcmdr("variable.list.height"),
                           selectmode=getRcmdr("multiple.select.mode"), export=FALSE)
    levelsScrollbar <- ttkscrollbar(levelsFrame, command=function(...) tkyview(levelsBox, ...))
    tkconfigure(levelsBox, yscrollcommand=function(...) tkset(levelsScrollbar, ...))
    tkselection.set(levelsBox, 0)

    onSelectTimeVar <- function() {
        var <- getSelection(timeVarBox)

        if(var == .gettext("Date")) {
            tkdelete(formatEntry, "0", "end")
            tkinsert(formatEntry, "end", "%Y-%m-%d")
        }
        else if (var == .gettext("Time")) {
            tkdelete(formatEntry, "0", "end")
            tkinsert(formatEntry, "end", "%Y-%m-%d %H:%M")
        }
    }

    tkbind(timeVarBox$listbox, "<<ListboxSelect>>", onSelectTimeVar)

    onSelectGroup <- function() {
        var <- c("", vars)[as.numeric(tkcurselection(varsBox))+1]
        tkdelete(levelsBox, "0", "end")

        if(var == "")
            return()

        levs <- if(is.factor(meta(corpus, var)[[1]])) levels(droplevels(meta(corpus, var)[[1]]))
                else sort(unique(meta(corpus, var)[[1]]))
        for(level in as.character(levs[!is.na(levs)])) tkinsert(levelsBox, "end", level)

        tkselection.set(levelsBox, 0, "end")
    }

    tkbind(varsBox, "<<ListboxSelect>>", onSelectGroup)

    tclTerms <- tclVar("")
    entryTerms <- ttkentry(top, width="30", textvariable=tclTerms)

    radioButtons(name="what",
                 buttons=c("term.lev", "lev.term", "absolute"),
                 labels=c(.gettext("Term prevalence in level (\"% Term/Level\")"),
                          .gettext("Distribution of occurrences among levels (\"% Level/Term\")"),
                          .gettext("Number of occurrences in level per time unit")),
                 title=.gettext("Measure:"),
                 right.buttons=FALSE)

    tclMean <- tclVar(0)
    meanButton <- tkcheckbutton(top, text=.gettext("Apply rolling mean"), variable=tclMean)

    tclWindow <- tclVar(7)
    spinWindow <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                             inc=1, textvariable=tclWindow,
                             validate="all", validatecommand=.validate.uint)

    tclTitle <- tclVar(.gettext("Temporal evolution of occurrences"))
    titleEntry <- ttkentry(top, width="30", textvariable=tclTitle)

    onCustom <- function() {
        timeVar <- getSelection(timeVarBox)
        groupVar <- c("", vars)[as.numeric(tkcurselection(varsBox))+1]
        termsList <- strsplit(tclvalue(tclTerms), " ")[[1]]
        what <- tclvalue(whatVariable)
        rollmean <- tclvalue(tclMean) == 1
        window <- as.numeric(tclvalue(tclWindow))
        title <- tclvalue(tclTitle)

        if(what == "lev.term" && nchar(groupVar) == 0) {
            .Message(message=.gettext("Plotting distribution of occurrences with only one curve does not make sense: all points would be 100%."),
                     type="error", parent=top)
            return()
        }

        if(length(termsList) > 1 && nchar(groupVar) > 0) {
            .Message(message=.gettext("Only one term can be used when a grouping variable is selected."),
                     type="error", parent=top)
            return()
        }

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
                            paste(wrongTerms, collapse=.gettext("\', \'"))), type="error", parent=top)
            return()
        }

        if(nchar(groupVar) > 0)
            groupLevs <- unique(meta(corpus, groupVar)[[1]])[as.numeric(tkcurselection(levelsBox))+1]

        format <- tclvalue(tclFormat)

        # Check that format is more or less correct before running the code
        time <- meta(corpus, timeVar)[[1]]
        time <- strptime(unique(time[!is.na(time)]), format)
        if(all(is.na(time))) {
            .Message(message=sprintf(.gettext("Incorrect time format or variable: no values of \"%s\" could be converted to a time index."), timeVar),
                     type="error", parent=top)
            return()
        }
        else if(any(is.na(time))) {
            .Message(message=sprintf(.gettext("Some values of \"%s\" could not be converted to a time index and will be missing."), timeVar),
                     type="warning", parent=top)
        }

        doItAndPrint(sprintf('time <- as.character(strptime(meta(corpus, "%s")[[1]], "%s"))', timeVar, format))

        if(nchar(groupVar) == 0) {
            doItAndPrint(sprintf('absTermFreqs <- as.table(rollup(dtm[, c("%s")], 1, time))', paste(termsList, collapse='", "')))
        }
        else {
            doItAndPrint(sprintf('absTermFreqs <- as.table(tapply(as.numeric(as.matrix(dtm[, c("%s")])), list(time, meta(corpus, "%s")[[1]]), sum))',
                                 termsList, groupVar))

            if(length(groupLevs) < length(unique(meta(corpus, var)[[1]])))
                doItAndPrint(sprintf('absTermFreqs <- absTermFreqs[,c(%s)]',
                                     paste(groupLevs, collapse='", "')))
        }

            doItAndPrint("names(dimnames(absTermFreqs)) <- NULL")

            # Compute %
            if(what == "term.lev") {
                doItAndPrint("termSeries <- zoo(absTermFreqs/c(tapply(row_sums(dtm), time, sum)), order.by=as.POSIXct(rownames(absTermFreqs)))")

                ylab <- .gettext("% of all terms")
            }
            else if (what == "lev.term") {
                doItAndPrint("termSeries <- zoo(prop.table(absTermFreqs, 2) * 100, order.by=as.POSIXct(rownames(absTermFreqs)))")
                ylab <- .gettext("% of occurrences")
            }
            else {
                doItAndPrint("termSeries <- zoo(absTermFreqs, order.by=as.POSIXct(rownames(absTermFreqs)))")
                ylab <- .gettext("Number of occurrences")
            }

        # We need to be sure we have a common time unit, i.e. we have a zooreg object
        if(length(termSeries) > 1 && !is.regular(termSeries))
            doItAndPrint('termSeries <- aggregate(termSeries, list(as.POSIXct(trunc(time, units(diff(time))))), regular=TRUE)')


        # For some reason, computing this after merging returns 24 "hours" instead of 1 "day" as unit
        unitsLoc <- c(secs=.gettext("per second"), mins=.gettext("per minute"), hours=.gettext("per hour"),
                      days=.gettext("per day"), weeks=.gettext("per week"))
        unit <- unitsLoc[units(diff(time(termSeries)))]

        # Trick to get a strictly regular time series with NA where no document was observed
        # difftime chooses the unit so that all differences are > 1, which is what we want
        if(length(termSeries) > 1 && !is.regular(termSeries, strict=TRUE)) {
            # seq() will specify different if we pass "day"
            byUnit <- units(diff(time(termSeries)))
            if(byUnit == "days") byUnit <- "DSTday"

            doItAndPrint(sprintf('termSeries <- merge(termSeries, zoo(, seq(start(termSeries), end(termSeries), "%s")), fill=%s)',
                                 byUnit, if(what == "absolute") "0" else "NaN"))
        }

        if(rollmean) {
            if(window >= NROW(termSeries))
                .Message(message=.gettext("Chosen roll mean window is longer than the range of the time variable, rolling mean was not applied."),
                         type="warning", parent=top)
            else
                # For percents, the days with no observation get 0/0 == NaN, and we need to skip them
                doItAndPrint(sprintf('termSeries <- rollapply(termSeries, %s, align="left", mean, na.rm=TRUE)', window))
        }

        # If only one time point is present, plotting always fails
        if(length(termSeries) > 1)
            doItAndPrint(sprintf('xyplot(termSeries, superpose=TRUE, xlab="", ylab="%s", main="%s", auto.key=%s, par.settings=simpleTheme(lwd=1.5))',
                                 paste(ylab, unit), title,
                                 if(NCOL(termSeries) > 1) 'TRUE' else "NULL"))
        else
            .Message(.gettext("Only one time point present, no plot can be drawn."), "error", parent=top)

        doItAndPrint("rm(absTermFreqs, time)")

        setLastTable("termSeries", paste(title, " (", ylab, " ", unit, ")", sep=""))

        doItAndPrint("termSeries")

        activateMenus()
        tkfocus(CommanderWindow())
    }

    # Shut up R CMD check WARNING
    onClose <- NULL
    .customCloseHelp(helpSubject="varTimeSeriesDlg", custom.button=.gettext("Draw plot"))

    tkgrid(getFrame(timeVarBox), sticky="ewns", padx=c(0, 6), pady=6, row=0, rowspan=3)
    tkgrid(labelRcmdr(top, text=.gettext("Time format:")), padx=c(6, 0), pady=c(6, 0), sticky="w", row=0, column=1)
    tkgrid(formatEntry, sticky="w", padx=c(6, 0), row=1, column=1)
    tkgrid(labelRcmdr(top, text=.gettext("%Y: year - %m: month - %d: day\n%H: hour - %M: minute\nClick the \"Help\" button for more codes.")),
           sticky="w", row=2, column=1, pady=6)
    tkgrid(.titleLabel(varsFrame, text=.gettext("Group by variable:")),
           sticky="w", pady=c(12, 0))
    tkgrid(varsBox, varsScrollbar, sticky="ewns", pady=6)
    tkgrid(varsFrame, levelsFrame, sticky="ewns", pady=6)
    tkgrid(labelRcmdr(levelsFrame, text=.gettext("Only plot levels:")), sticky="w", pady=c(12, 0))
    tkgrid(levelsBox, levelsScrollbar, sticky="ewns", pady=6)
    tkgrid(.titleLabel(top, text=.gettext("Terms to show (space-separated):")),
           sticky="w", pady=c(6, 0))
    tkgrid(entryTerms, sticky="w", pady=6, columnspan=2)
    tkgrid(whatFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(.titleLabel(top, text=.gettext("Rolling mean:")),
           sticky="w", pady=c(6, 0))
    tkgrid(meanButton, sticky="w")
    tkgrid(labelRcmdr(top, text=.gettext("Time window for mean (in time units):")), spinWindow, sticky="w",
           padx=6, pady=c(0, 6))
    tkgrid(.titleLabel(top, text=.gettext("Title:")),
           sticky="w", pady=c(6, 0))
    tkgrid(titleEntry, sticky="w", padx=6, pady=c(0, 6), columnspan=2)
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)
    dialogSuffix(focus=timeVarBox$listbox, onOK=onCustom, onCancel=onClose)
}
