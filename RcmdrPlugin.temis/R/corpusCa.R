runCorpusCa <- function(corpus, dtm=NULL, variables=NULL, sparsity=0.9, ...) {
    if(is.null(dtm))
        dtm<-DocumentTermMatrix(corpus)

    if(!all(variables %in% colnames(meta(corpus))))
        stop("All items of 'variables' should be meta-data variables of the corpus.")

    # Save old meta-data now to check what is lost when skipping documents
    oldMeta<-meta<-meta(corpus)[colnames(meta(corpus)) != "MetaID"]

    # removeSparseTerms() does not accept 1
    if(sparsity < 1)
        dtm<-removeSparseTerms(dtm, sparsity)

    invalid<-which(apply(dtm,1,sum)==0)
    if(length(invalid) > 0) {
        dtm<-dtm[-invalid, , drop=FALSE]
        meta<-oldMeta[-invalid, , drop=FALSE]
        corpus<-corpus[-invalid]
    }

    ndocs<-nrow(dtm)
    nterms<-ncol(dtm)

    if(ndocs <= 1 || nterms <= 1) {
        .Message(.gettext("Please increase the value of the 'sparsity' parameter so that at least two documents and two terms are retained."),
                 type="error")
        return()
    }

    if(length(invalid) > 0) {
        Message(paste(.gettext("Documents skipped from correspondence analysis:\n"),
                      paste(names(invalid), collapse=", ")),
                type="note")

        .Message(sprintf(.gettext("%i documents have been skipped because they do not include any occurrence of the terms retained in the final document-term matrix. Their list is available in the \"Messages\" area.\n\nIncrease the value of the 'sparsity' parameter if you want to include them."),
                         length(invalid)),
                 type="info")
    }

    skippedVars<-character()
    skippedLevs<-character()
    origVars<-character()

    dupLevels<-any(duplicated(unlist(lapply(meta, function(x) substr(unique(as.character(x[!is.na(x)])), 0, 30)))))


    varDtm <- NULL

    # Create mean dummy variables as rows
    # Keep in sync with showCorpusClustering()

    # Just in case variables have common levels, and are truncated to the same string
    vars <- colnames(meta)
    vars10<-make.unique(substr(vars, 0, 10))
    vars20<-make.unique(substr(vars, 0, 20))

    if(ncol(meta) > 0) {
        for(i in 1:ncol(meta)) {
            var<-vars[i]
            levs<-levels(factor(meta[,i]))
            totNLevels<-nlevels(factor(oldMeta[,i]))

            if(length(levs) == 0) {
                skippedVars<-c(skippedVars, var)
                next
            }
            else if(length(levs) < totNLevels) {
                skippedLevs<-c(skippedLevs, var)
            }

            # suppressWarnings() is used because rollup() warns when NAs are present
            suppressWarnings(mat<-rollup(dtm[1:ndocs, , drop=FALSE], 1, meta[i]))

            # If only one level is present, don't add the level name (e.g. YES),
            # except if all values are the same (in which case variable is useless but is more obvious that way)
            if(totNLevels == 1 && any(is.na(meta[,i])))
                rownames(mat)<-vars20[i]
            # In case of ambiguous levels of only numbers in levels, add variable names everywhere
            else if(dupLevels || !any(is.na(suppressWarnings(as.numeric(levs)))))
                rownames(mat)<-make.unique(paste(vars10[i], substr(levs, 0, 30)))
            else # Most general case: no need to waste space with variable names
                rownames(mat)<-substr(levs, 0, 30)

            varDtm<-rbind(varDtm, mat)
            origVars<-c(origVars, rep(var, nrow(mat)))
        }
    }

    if(!is.null(variables) && sum(origVars %in% variables) < 2) {
        .Message(.gettext("Please select active variables so that at least two levels are present in the retained documents."),
                type="error")
        return()
    }

    Message(sprintf(.gettext("Running correspondence analysis using %i documents, %i terms and %i variables."),
                    ndocs, nterms, ncol(meta)),
            type="note")

    if(length(skippedVars) > 0)
        msg1 <- sprintf(.gettext("Variable(s) %s have been skipped since it contains only missing values for retained documents."),
                        paste(skippedVars, collapse=", "))
    else
        msg1 <- ""

    if(length(skippedLevs) > 0)
        msg2 <- sprintf(.gettext("Some levels of variable(s) %s have been skipped since they contain only missing values for retained documents."),
                        paste(skippedLevs, collapse=", "))
    else
        msg2 <- ""

    if(length(skippedVars) > 0 && length(skippedLevs) > 0)
        .Message(paste(msg1, "\n\n", msg2), "info")
    else if(length(skippedVars) > 0)
        .Message(msg1, "info")
    else if(length(skippedLevs) > 0)
        .Message(msg2, "info")

    newDtm <- as.matrix(rbind(dtm, varDtm))

    if(!is.null(variables))
        obj <- ca(newDtm, suprow=c(1:nrow(dtm), nrow(dtm) + which(!origVars %in% variables)), ...)
    else if(nrow(newDtm) - ndocs > 0)
        obj <- ca(newDtm, suprow=(ndocs+1):nrow(newDtm), ...)
    else
        obj <- ca(newDtm, ...)

    if(nrow(newDtm) - ndocs > 0) {
        obj$rowvars <- seq.int(ndocs + 1, nrow(newDtm))
        names(obj$rowvars) <- origVars
    }

    # This is used by corpusClustDlg() when computing distances between documents using dist()
    rownames(obj$rowcoord) <- obj$rownames
    rownames(obj$colcoord) <- obj$colnames

    attr(obj, "sparsity") <- sparsity

    obj
}

corpusCaDlg <- function() {
    initializeDialog(title=.gettext("Run Correspondence Analysis"))

    labelNDocs <- labelRcmdr(top)

    labels <- c(.gettext("(Terms present in at least %s documents will be retained in the analysis.)"),
                .gettext("(All terms will be retained in the analysis.)"))

    tkconfigure(labelNDocs, width=max(nchar(labels)))

    updateNDocs <- function(value) {
        ndocs <- ceiling((1 - as.numeric(value)/100) * nrow(dtm))

        if(ndocs > 1)
            tkconfigure(labelNDocs, text=sprintf(labels[1], ndocs))
        else
            tkconfigure(labelNDocs, text=labels[2])
    }

    vars <- c(.gettext("None (run analysis on full matrix)"), colnames(meta(corpus)))
    varBox <- variableListBox(top, vars,
                              selectmode="multiple",
                              title=.gettext("Aggregate document-term matrix by variables:"),
                              initialSelection=0)

    tclSparsity <- tclVar(100 - ceiling(1/nrow(dtm) * 100))
    spinSparsity <- tkwidget(top, type="spinbox", from=0, to=100,
                             inc=0.1, textvariable=tclSparsity,
                             validate="all", validatecommand=function(P) .validate.unum(P, fun=updateNDocs))
    updateNDocs(tclvalue(tclSparsity))


    tclDim <- tclVar(5)
    sliderDim <- tkscale(top, from=1, to=30,
                         showvalue=TRUE, variable=tclDim,
	                     resolution=1, orient="horizontal")


    onOK <- function() {
        sparsity <- as.numeric(tclvalue(tclSparsity))
        vars <- getSelection(varBox)
        dim <- as.numeric(tclvalue(tclDim))

        if(is.na(sparsity) || sparsity <= 0 || sparsity > 100) {
            .Message(.gettext("Please specify a sparsity value between 0 (excluded) and 100%."), type="error")
            return()
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        if(ncol(meta(corpus)[colnames(meta(corpus)) != "MetaID"]) == 0)
            Message(message=.gettext("No corpus variables have been set. Use Text mining->Manage corpus->Set corpus variables to add them."),
                    type="note")

        if(length(vars) == 1 && vars[1] == .gettext("None (run analysis on full matrix)"))
            doItAndPrint(sprintf("corpusCa <- runCorpusCa(corpus, dtm, sparsity=%s, nd=%i)", sparsity/100, dim))
        else
            doItAndPrint(sprintf('corpusCa <- runCorpusCa(corpus, dtm, c("%s"), sparsity=%s, nd=%i)',
                                  paste(vars, collapse='", "'), sparsity/100, dim))

        if(!is.null(corpusCa)) {
            setLastTable("corpusCa", .gettext("Correspondence analysis"))

            showCorpusCaDlg()
        }

        activateMenus()

        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject=corpusCaDlg)
    tkgrid(getFrame(varBox), columnspan=2, sticky="we", pady=6)
    tkgrid(labelRcmdr(top, text=.gettext("Remove terms missing from more than (% of documents):")),
           spinSparsity, sticky="sew", pady=6)
    tkgrid(labelNDocs, sticky="sw", pady=6, columnspan=2)
    tkgrid(labelRcmdr(top, text=.gettext("Number of dimensions to retain:")),
           sliderDim, sticky="sew", pady=6)
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix()
}

