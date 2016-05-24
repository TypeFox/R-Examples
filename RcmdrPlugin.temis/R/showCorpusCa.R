## Helper functions to plot a CA restricted to most contributive points

# Get rows contributions to the dimensions of a CA
rowCtr <- function(obj, dim) {
    I <- dim(obj$rowcoord)[1]
    # nd can be lower than actual number of dimensions if matrix does not require them
    K <- min(obj$nd, ncol(obj$rowcoord))
    svF <- matrix(rep(obj$sv[1:K], I), I, K, byrow=TRUE)
    rpc <- obj$rowcoord[,1:K] * svF
    obj$rowmass * rowSums(rpc[,dim, drop=FALSE]^2) / sum(obj$sv[dim]^2)
}

# Get columns contributions to the dimensions of a CA
colCtr <- function(obj, dim) {
    J <- dim(obj$colcoord)[1]
    # nd can be lower than actual number of dimensions if matrix does not require them
    K <- min(obj$nd, ncol(obj$colcoord))
    svG <- matrix(rep(obj$sv[1:K], J), J, K, byrow=TRUE)
    cpc <- obj$colcoord[,1:K] * svG
    obj$colmass * rowSums(cpc[,dim, drop=FALSE]^2) / sum(obj$sv[dim]^2)
}

# Restrain CA to a subset of rows (only for plotting!)
rowSubsetCa <- function(obj, indices) {
    ret <- obj
    for (i in 4:7) {
        ret[[i]] <- list()
        ret[[i]] <- obj[[i]][indices]
    }
    del <- which(!(1:nrow(obj$rowcoord) %in% indices))
    ret$rowsup <- as.numeric(sapply(obj$rowsup[obj$rowsup %in% indices],
                                    function(x) x - sum(del < x)))
    ret$rowcoord <- matrix()
    ret$rowcoord <- obj$rowcoord[indices,,drop=FALSE]
    ret$rownames <- obj$rownames[indices]
    ret
}

# Restrain CA to a subset of columns (only for plotting!)
colSubsetCa <- function(obj, indices) {
    ret <- obj
    for (i in 9:12) {
        ret[[i]] <- list()
        ret[[i]] <- obj[[i]][indices]
    }
    ret$colsup <- ret$colsup[ret$colsup %in% indices]
    ret$colsup <- as.numeric(lapply(obj$colsup, function(x) x - sum(indices < x)))
    ret$colcoord <- matrix()
    ret$colcoord <- obj$colcoord[indices,,drop=FALSE]
    ret$colnames <- obj$colnames[indices]
    ret
}

showCorpusCa <- function(corpusCa, dim=1, ndocs=10, nterms=10) {
    objects <- .getCorpusWindow()
    window <- objects$window
    txt <- objects$txt
    listbox <- objects$listbox

    tkwm.title(window, .gettext("Correspondence Analysis"))

    # Are active rows documents, rather than variables?
    actDocs<-length(corpusCa$rowvars) == length(corpusCa$rowsup)

    if(actDocs) {
        tndocs <- nrow(corpusCa$rowcoord) - length(corpusCa$rowsup)
        tnterms <- nrow(corpusCa$colcoord)
        tnvars <- length(unique(names(corpusCa$rowvars)))

        tkinsert(txt, "end", sprintf(.gettext("Correspondence analysis of %i documents, %i terms and %i supplementary variables."),
                                     tndocs, tnterms, tnvars),
                 "body")
    }
    else {
        tnactvars <- length(unique(names(corpusCa$rowvars)[!corpusCa$rowvars %in% corpusCa$rowsup]))
        tndocs <-  nrow(corpusCa$rowcoord) - length(corpusCa$rowvars)
        tnterms <- nrow(corpusCa$colcoord)
        tnsupvars <- length(unique(names(corpusCa$rowvars)[corpusCa$rowvars %in% corpusCa$rowsup]))

        tkinsert(txt, "end", sprintf(.gettext("Correspondence analysis of %i active variable(s) (aggregating %i documents), %i terms and %i supplementary variable(s)."),
                                     tnactvars, tndocs, tnterms, tnsupvars),
                 "body")
    }

    tkinsert(txt, "end", "\n\n", "heading")

    mark <- 0

    titles <- c(.gettext("Position"), .gettext("Contribution (%)"), .gettext("Quality (%)"))

    tkinsert(txt, "end", paste(.gettext("Axes summary:"), "\n", sep=""), "heading")
    tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
    tkinsert(listbox, "end", .gettext("Axes summary"))
    mark <- mark + 1

    nd <- min(length(corpusCa$sv), corpusCa$nd)
    values <- 100 * (corpusCa$sv[1:nd]^2)/sum(corpusCa$sv^2)
    values2 <- cumsum(values)
    val <- rbind(values, values2)
    rownames(val) <- c(.gettext("Inertia (%)"), .gettext("Cumulated inertia (%)"))
    colnames(val) <- seq.int(ncol(val))
    names(dimnames(val)) <- c("", .gettext("Axis"))
    tkinsert(txt, "end", paste(capture.output(val), collapse="\n"), "fixed")

    # Contributions to both axes
    if(actDocs) {
        rows <- order(rowCtr(corpusCa, dim), decreasing=TRUE)[1:ndocs]
        rows <- rows[!rows %in% corpusCa$rowsup]
    }
    else {
        rows <- corpusCa$rowvars[!corpusCa$rowvars %in% corpusCa$rowsup]
    }

    cols <- order(colCtr(corpusCa, dim), decreasing=TRUE)[1:nterms]
    cols <- cols[!cols %in% corpusCa$colsup]


    for(j in 1:length(dim)) {
        # Per-axis contributions, decreasing order
        rowsCtr <- rowCtr(corpusCa, dim[j])
        rows <- rows[order(rowsCtr[rows], decreasing=TRUE)]

        colsCtr <- colCtr(corpusCa, dim[j])
        cols <- cols[order(colsCtr[cols], decreasing=TRUE)]

        tkinsert(txt, "end",
                 paste("\n\n", sprintf(.gettext("Most contributive terms on negative side of axis %i:"), dim[j]), "\n", sep=""),
                 "heading")
        tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
        tkinsert(listbox, "end", sprintf(.gettext("Axis %i - Negative Side:"), dim[j]))
        tkitemconfigure(listbox, mark, background="grey")
        mark <- mark + 1

        negcols <- cols[corpusCa$colcoord[cols, dim[j]] < 0]
        if(length(negcols) == 0) {
            tkinsert(txt, "end",
                     sprintf(.gettext("None among the %i most contributive terms."), nterms), "body")
        }
        else {
            df <- data.frame(row.names=corpusCa$colnames[negcols],
                             corpusCa$colcoord[negcols, dim[j]] * corpusCa$sv[dim[j]],
                             colsCtr[negcols] * 100,
                             (corpusCa$colcoord[negcols, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$coldist[negcols])^2 * 100)
            colnames(df) <- titles

            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")
        }


        if(!actDocs) {
            negrows <- rows[corpusCa$rowcoord[rows, dim[j]] <= 0]

            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Active levels on negative side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")
            df <- data.frame(row.names=corpusCa$rownames[negrows],
                             corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]],
                             rowsCtr[negrows] * 100,
                             (corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[negrows])^2 * 100)
            colnames(df) <- titles

            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")


            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Most extreme documents on negative side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")

            int <- setdiff(1:nrow(corpusCa$rowcoord), corpusCa$rowvars)
            negrows <- int[corpusCa$rowcoord[int, dim[j]] < 0]
            negrows <- negrows[order(corpusCa$rowcoord[negrows, dim[j]])[1:ndocs]]
            negrows <- negrows[!is.na(negrows)]

            df <- data.frame(row.names=corpusCa$rownames[negrows],
                             corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]],
                             (corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[negrows])^2 * 100)
            colnames(df) <- titles[-2]
        }
        else {
            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Most contributive documents on negative side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")

            negrows <- rows[corpusCa$rowcoord[rows, dim[j]] < 0]

            df <- data.frame(row.names=corpusCa$rownames[negrows],
                             corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]],
                             rowsCtr[negrows] * 100,
                             (corpusCa$rowcoord[negrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[negrows])^2 * 100)
            colnames(df) <- titles
        }

        if(length(negrows) == 0) {
            tkinsert(txt, "end",
                     sprintf(.gettext("None among the %i most contributive documents."), ndocs), "body")
        }
        else {
            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")

            tkinsert(txt, "end", "\n", "heading")

            for(i in negrows) {
                # We need to use IDs rather than indexes to access documents in the corpus
                # since some documents may have been skipped in the CA
                id <- corpusCa$rownames[i]
                tkinsert(txt, "end", paste("\n", id, "\n", sep=""),
                         "articlehead")
                tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
                mark <- mark + 1
                tkinsert(listbox, "end", id)

                origin <- meta(corpus[[id]], "origin")
                date <- meta(corpus[[id]], "datetimestamp")
                if(length(origin) > 0 && length(date) > 0)
                    tkinsert(txt, "end", paste(origin, " - ", date, "\n", sep=""), "details")
                else if(length(origin) > 0)
                    tkinsert(txt, "end", paste(origin, "\n", sep=""), "details")
                else if(length(origin) > 0)
                    tkinsert(txt, "end", paste(date, "\n", sep=""), "details")

                if(length(origin) > 0 || length(date) > 0)
                    tkinsert(txt, "end", "\n", "small")

                tkinsert(txt, "end", paste(paste(corpus[[id]], collapse="\n"), "\n"), "body")
            }
        }


        tkinsert(txt, "end",
                 paste("\n\n", sprintf(.gettext("Most contributive terms on positive side of axis %i:"),
                                       dim[j]), "\n", sep=""),
                 "heading")
        tkinsert(listbox, "end", sprintf(.gettext("Axis %i - Positive Side:"), dim[j]))
        tkitemconfigure(listbox, mark, background="grey")
        tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
        mark <- mark + 1

        poscols <- cols[corpusCa$colcoord[cols, dim[j]] >= 0]
        if(length(poscols) == 0) {
            tkinsert(txt, "end",
                     sprintf(.gettext("None among the %i most contributive terms."), nterms), "body")
        }
        else {
            df <- data.frame(row.names=corpusCa$colnames[poscols],
                             corpusCa$colcoord[poscols, dim[j]] * corpusCa$sv[dim[j]],
                             colsCtr[poscols] * 100,
                             (corpusCa$colcoord[poscols, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$coldist[poscols])^2 * 100)
            colnames(df) <- titles

            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")
        }

        if(!actDocs) {
            posrows <- rows[corpusCa$rowcoord[rows, dim[j]] > 0]

            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Active levels on positive side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")
            df <- data.frame(row.names=corpusCa$rownames[posrows],
                             corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]],
                             rowsCtr[posrows] * 100,
                             (corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[posrows])^2 * 100)
            colnames(df) <- titles

            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")

            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Most extreme documents on positive side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")

            int <- setdiff(1:nrow(corpusCa$rowcoord), corpusCa$rowvars)
            posrows <- int[corpusCa$rowcoord[int, dim[j]] > 0]
            posrows <- posrows[order(corpusCa$rowcoord[posrows, dim[j]], decreasing=TRUE)[1:ndocs]]
            posrows <- posrows[!is.na(posrows)]

            df <- data.frame(row.names=corpusCa$rownames[posrows],
                             corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]],
                             (corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[posrows])^2 * 100)
            colnames(df) <- titles[-2]
        }
        else {
            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Most contributive documents on positive side of axis %i:"), dim[j]), "\n", sep=""),
                     "heading")

            posrows <- rows[corpusCa$rowcoord[rows, dim[j]] > 0]

            df <- data.frame(row.names=corpusCa$rownames[posrows],
                             corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]],
                             rowsCtr[posrows] * 100,
                             (corpusCa$rowcoord[posrows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[posrows])^2 * 100)
            colnames(df) <- titles
        }

        if(length(posrows) == 0) {
            tkinsert(txt, "end",
                     sprintf(.gettext("None among the %i most contributive documents."), ndocs), "body")
        }
        else {
            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")

            tkinsert(txt, "end", "\n", "heading")

            for(i in posrows) {
                # We need to use IDs rather than indexes to access documents in the corpus
                # since some documents may have been skipped in the CA
                id <- corpusCa$rownames[i]
                tkinsert(txt, "end", paste("\n", id, "\n", sep=""),
                         "articlehead")
                tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
                mark <- mark + 1
                tkinsert(listbox, "end", id)

                origin <- meta(corpus[[id]], "origin")
                date <- meta(corpus[[id]], "datetimestamp")
                if(length(origin) > 0 && length(date) > 0)
                    tkinsert(txt, "end", paste(origin, " - ", date, "\n", sep=""), "details")
                else if(length(origin) > 0)
                    tkinsert(txt, "end", paste(origin, "\n", sep=""), "details")
                else if(length(origin) > 0)
                    tkinsert(txt, "end", paste(date, "\n", sep=""), "details")

                if(length(origin) > 0 || length(date) > 0)
                    tkinsert(txt, "end", "\n", "small")

                tkinsert(txt, "end", paste(paste(corpus[[id]], collapse="\n"), "\n"), "body")
            }
        }


        # Supplementary variables start after supplementary documents, if any
        suprows <- intersect(corpusCa$rowvars, corpusCa$rowsup)
        if(length(suprows) > 0) {
            tkinsert(txt, "end",
                      paste("\n\n", sprintf(.gettext("Situation of supplementary variables on axis %i:"), dim[j]), "\n", sep=""),
                     "heading")
            tkinsert(listbox, "end", sprintf(.gettext("Axis %i - Variables"), dim[j]))
            tkitemconfigure(listbox, mark, background="grey")
            tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
            mark <- mark + 1

            df <- data.frame(row.names=corpusCa$rownames[suprows],
                             corpusCa$rowcoord[suprows, dim[j]] * corpusCa$sv[dim[j]],
                             (corpusCa$rowcoord[suprows, dim[j]] * corpusCa$sv[dim[j]] / corpusCa$rowdist[suprows])^2 * 100)
            colnames(df) <- titles[-2]

            tkinsert(txt, "end", paste(capture.output(format(df)), collapse="\n"), "fixed")
        }
    }

    tkraise(window)
}

showCorpusCaDlg <- function() {
    if(!exists("corpusCa") || !class(corpusCa) == "ca") {
        .Message(message=.gettext("Please run a correspondence analysis on the corpus first."),
                 type="error")
        return()
    }

    initializeDialog(title=.gettext("Show Correspondence Analysis"))

    actDocs<-length(corpusCa$rowvars) == length(corpusCa$rowsup)

    dimFrame <- tkframe(top)
    tkgrid(.titleLabel(dimFrame, text=.gettext("Dimensions to plot:")),
           sticky="s")
    tclXDim <- tclVar(1)
    tclYDim <- tclVar(2)
    nd <- min(length(corpusCa$sv), corpusCa$nd)
    xSlider <- tkscale(dimFrame, from=1, to=nd,
                      showvalue=TRUE, variable=tclXDim,
		      resolution=1, orient="horizontal")
    ySlider <- tkscale(dimFrame, from=1, to=nd,
                      showvalue=TRUE, variable=tclYDim,
		      resolution=1, orient="horizontal")

    if(!actDocs) {
        checkBoxes(frame="labelsFrame",
                   boxes=c("varLabels", "termLabels"),
                   initialValues=c(1, 1),
                   labels=c(.gettext("Variables"), .gettext("Terms")),
                   title=.gettext("Draw labels for:"))

        checkBoxes(frame="pointsFrame",
                   boxes=c("varPoints", "termPoints"),
                   initialValues=c(1, 1),
                   labels=c(.gettext("Variables"), .gettext("Terms")),
                   title=.gettext("Draw point symbols for:"))
    }
    else if(length(corpusCa$rowsup) == 0) {
        checkBoxes(frame="labelsFrame",
                   boxes=c("docLabels", "termLabels"),
                   initialValues=c(0, 1),
                   labels=c(.gettext("Documents"), .gettext("Terms")),
                   title=.gettext("Draw labels for:"))

        checkBoxes(frame="pointsFrame",
                   boxes=c("docPoints", "termPoints"),
                   initialValues=c(0, 1),
                   labels=c(.gettext("Documents"), .gettext("Terms")),
                   title=.gettext("Draw point symbols for:"))
    }
    else {
        checkBoxes(frame="labelsFrame",
                   boxes=c("varLabels", "docLabels", "termLabels"),
                   initialValues=c(0, 0, 1),
                   labels=c(.gettext("Variables"), .gettext("Documents"), .gettext("Terms")),
                   title=.gettext("Draw labels for:"))

        checkBoxes(frame="pointsFrame",
                   boxes=c("varPoints", "docPoints", "termPoints"),
                   initialValues=c(0, 0, 1),
                   labels=c(.gettext("Variables"), .gettext("Documents"), .gettext("Terms")),
                   title=.gettext("Draw point symbols for:"))
    }

    vars <- colnames(meta(corpus))

    if(actDocs)
        selection <- (1:length(vars)) - 1
    else
        selection <- match(unique(names(corpusCa$rowvars[!corpusCa$rowvars %in% corpusCa$rowsup])), vars) - 1

    varBox <- variableListBox(top, vars,
                              selectmode="multiple",
                              title=.gettext("Variables to plot:"),
                              initialSelection=selection)
    tkbind(varBox$listbox, "<<ListboxSelect>>", function(...) tclvalue(varLabelsVariable) <- 1)

    nFrame <- tkframe(top)
    tkgrid(.titleLabel(nFrame, text=.gettext("Number of items to plot:")),
           sticky="s")
    tclNDocs <- tclVar(25)
    tclNTerms <- tclVar(25)
    spinDocs <- tkwidget(top, type="spinbox", from=1, to=nrow(corpusCa$rowcoord)-length(corpusCa$rowsup),
                         inc=1, textvariable=tclNDocs,
                         validate="all", validatecommand=.validate.uint)

    spinTerms <- tkwidget(top, type="spinbox", from=1, to=nrow(corpusCa$colcoord)-length(corpusCa$colsup),
                          inc=1, textvariable=tclNTerms,
                          validate="all", validatecommand=.validate.uint)

    ctrDimVariable <- tclVar("xyDim")
    ctrDimFrame <- tkframe(top)
    ctrDim1 <- ttkradiobutton(ctrDimFrame, variable=ctrDimVariable, value="xyDim", text=.gettext("Both axes"))
    ctrDim2 <- ttkradiobutton(ctrDimFrame, variable=ctrDimVariable, value="xDim", text=.gettext("Horizontal axis"))
    ctrDim3 <- ttkradiobutton(ctrDimFrame, variable=ctrDimVariable, value="yDim", text=.gettext("Vertical axis"))


    onCustom <- function() {
        x <- tclvalue(tclXDim)
        y <- tclvalue(tclYDim)
        docLabels <- if(!actDocs) FALSE else tclvalue(docLabelsVariable) == 1
        termLabels <- tclvalue(termLabelsVariable) == 1
        varLabels <- if(actDocs && length(corpusCa$rowsup) == 0) FALSE else tclvalue(varLabelsVariable) == 1
        docPoints <- if(!actDocs) FALSE else tclvalue(docPointsVariable) == 1
        termPoints <- tclvalue(termPointsVariable) == 1
        varPoints <- if(actDocs && length(corpusCa$rowsup) == 0) FALSE else tclvalue(varPointsVariable) == 1
        vars <- getSelection(varBox)
        nDocs <- as.integer(tclvalue(tclNDocs))
        nTerms <- as.integer(tclvalue(tclNTerms))
        ctrDim <- switch(tclvalue(ctrDimVariable), xyDim=paste("c(", x, ", ", y, ")", sep=""), xDim=x, yDim=y)

        if(nDocs > nrow(corpusCa$rowcoord) - length(corpusCa$rowsup))
            nDocs <- nrow(corpusCa$rowcoord) - length(corpusCa$rowsup)

        if(nTerms > nrow(corpusCa$colcoord) - length(corpusCa$colsup))
            nTerms <- nrow(corpusCa$colcoord) - length(corpusCa$colsup)

        if(!(docLabels || termLabels || varLabels || docPoints || termPoints || varPoints)) {
            .Message(.gettext("Please select something to plot."), "error", parent=top)

            return()
        }

        setBusyCursor()
        on.exit(setIdleCursor())

        doItAndPrint(sprintf("showCorpusCa(corpusCa, %s, %s, %s)", ctrDim, nDocs, nTerms))


        if(actDocs) {
            if((docLabels || docPoints) && (varLabels || varPoints)) {
                rowWhat <- "all"
                varIndexes <- corpusCa$rowvars[names(corpusCa$rowvars) %in% vars]
                doItAndPrint(paste("plottingCa <- rowSubsetCa(corpusCa, c(order(rowCtr(corpusCa, ", ctrDim,
                                   "), decreasing=TRUE)[1:", nDocs, "], ", paste(varIndexes, collapse=", "), "))", sep=""))
            }
            else if(docLabels || docPoints) {
                rowWhat <- "active"
                doItAndPrint(paste("plottingCa <- rowSubsetCa(corpusCa, order(rowCtr(corpusCa, ", ctrDim,
                                   "), decreasing=TRUE)[1:", nDocs, "])", sep=""))
            }
            else if(varLabels || varPoints) {
                rowWhat <- "passive"
                varIndexes <- corpusCa$rowvars[names(corpusCa$rowvars) %in% vars]
                doItAndPrint(paste("plottingCa <- rowSubsetCa(corpusCa, c(",
                                   paste(varIndexes, collapse=", "), "))", sep=""))
            }
            else {
                rowWhat <- "none"
            }

            if(((docPoints || docLabels) && (varPoints || varLabels)) && docLabels != varLabels)
                Message(.gettext("Plotting documents and variables at the same time currently forces labels to be drawn for both or none."),
                         "note")

            rowActivePoints<-if(docPoints) 16 else NA
            rowSupPoints<-if(varPoints) 1 else NA
            rowActiveColor<-"black"
            rowSupColor<-"blue"
            rowActiveFont<-3
            rowSupFont<-4
        }
        else {
            if(varLabels || varPoints) {
                rowWhat <- "all"
                varIndexes <- corpusCa$rowvars[names(corpusCa$rowvars) %in% vars]
                doItAndPrint(paste("plottingCa <- rowSubsetCa(corpusCa, c(",
                                   paste(varIndexes, collapse=", "), "))", sep=""))
            }
            else {
                rowWhat <- "none"
            }

            rowActivePoints<-if(varPoints) 1 else NA
            rowSupPoints<-if(varPoints) 16 else NA
            rowActiveColor<-rowSupColor<-"blue"
            rowActiveFont<-rowSupFont<-4
        }

        if(termLabels || termPoints) {
            colWhat <- "all"
            if(docLabels || docPoints || varLabels || varPoints)
                doItAndPrint(paste("plottingCa <- colSubsetCa(plottingCa, order(colCtr(corpusCa, ", ctrDim,
                                   "), decreasing=TRUE)[1:", nTerms, "])", sep=""))
            else
                doItAndPrint(paste("plottingCa <- colSubsetCa(corpusCa, order(colCtr(corpusCa, ", ctrDim,
                                   "), decreasing=TRUE)[1:", nTerms, "])", sep=""))
        }
        else {
            colWhat <- "none"
        }

        doItAndPrint(sprintf('plotCorpusCa(plottingCa, dim=c(%s, %s), what=c("%s", "%s"), labels=c(%i, %i), pch=c(%s, %s, %s, NA), col.text=c("%s", "%s", "black", "red"), font=c(%i, %i, 1, 2), mass=TRUE, xlab="%s", ylab="%s")',
                              x, y, rowWhat, colWhat,
                              if(docLabels  || varLabels) 2 else 0, if(termLabels) 2 else 0,
                              rowActivePoints, rowSupPoints, if(termPoints) 17 else NA,
                              rowActiveColor, rowSupColor, rowActiveFont, rowSupFont,
                              sprintf(.gettext("Dimension %s (%.1f%%)"), x, 100 * corpusCa$sv[as.integer(x)]^2/sum(corpusCa$sv^2)),
                              sprintf(.gettext("Dimension %s (%.1f%%)"), y, 100 * corpusCa$sv[as.integer(y)]^2/sum(corpusCa$sv^2))))


        activateMenus()
    }

    # Shut up R CMD check WARNING
    onClose <- NULL
    .customCloseHelp(helpSubject="showCorpusCaDlg", custom.button=.gettext("Show"))

    tkgrid(labelRcmdr(dimFrame, text=.gettext("Horizontal axis:")), xSlider, sticky="w")
    tkgrid(labelRcmdr(dimFrame, text=.gettext("Vertical axis:")), ySlider, sticky="w")
    tkgrid(dimFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(labelsFrame, pointsFrame, sticky="w", pady=6, padx=c(0, 6))
    if(!actDocs || length(corpusCa$rowsup) > 0)
        tkgrid(getFrame(varBox), columnspan=2, sticky="we", pady=6)
    if(actDocs)
        tkgrid(labelRcmdr(nFrame, text=.gettext("Documents:")), spinDocs, sticky="w")
    tkgrid(labelRcmdr(nFrame, text=.gettext("Terms:")), spinTerms, sticky="w")
    tkgrid(nFrame, sticky="w", pady=6, columnspan=2)
    tkgrid(labelRcmdr(ctrDimFrame, text=.gettext("Most contributive to:")), sticky="w", columnspan=2, pady=6)
    tkgrid(ctrDim1, ctrDim2, ctrDim3, sticky="w", pady=6)
    tkgrid(ctrDimFrame, sticky="w", pady=6, columnspan=2)
    tkgrid.columnconfigure(ctrDimFrame, "all", uniform="a")
    tkgrid(buttonsFrame, sticky="ew", pady=6, columnspan=2)

    # We don't use dialogSuffix() itself because the dialog should not close,
    # and yet not call tkwait.window(): the plot does not draw correctly on Mac OS if we do.
    # The grab has also been removed since it prevents the user from scrolling the CA text window.
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onCustom)
    tkbind(top, "<Escape>", onClose)
    if (getRcmdr("double.click") && (!preventDoubleClick)) tkbind(window, "<Double-ButtonPress-1>", onCustom)
    tkwm.deiconify(top)
    tkfocus(top)
    if (getRcmdr("crisp.dialogs")) tclServiceMode(on=TRUE)
}

