showCorpusClustering <- function(corpusSubClust, ndocs=10, nterms=20, p=0.1, min.occ=5) {
    setBusyCursor()
    on.exit(setIdleCursor())

    objects <- .getCorpusWindow()
    window <- objects$window
    txt <- objects$txt
    listbox <- objects$listbox

    tkwm.title(window, .gettext("Hierarchical Clustering"))

    tnterms <- attr(corpusSubClust, "nterms")
    tndocs <- sum(sapply(corpusSubClust$lower, attr, "members"))

    if(is.null(attr(corpusSubClust, "caDim"))) {
        tkinsert(txt, "end", sprintf(.gettext("Hierarchical clustering of %i documents using %i terms (Ward's method with Chi-squared distance)."),
                                     tndocs, tnterms),
                 "body")
    }
    else {
        tndim <- attr(corpusSubClust, "caDim")

        tkinsert(txt, "end", sprintf(.gettext("Hierarchical clustering of %i documents using %i terms (Ward's method with Chi-squared distance based on the %i first dimensions of the correspondence analysis)."),
                                     tndocs, tnterms, tndim),
                 "body")
    }

    tkinsert(txt, "end", "\n\n", "heading")

    mark <- 0

    tkinsert(txt, "end", paste(.gettext("Clusters summary:"), "\n", sep=""), "heading")
    tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
    tkinsert(listbox, "end", .gettext("Clusters summary"))
    mark <- mark + 1

    val <- rbind(sapply(corpusSubClust$lower, attr, "members"),
                 sapply(corpusSubClust$lower, attr, "members")/sum(!is.na(meta(corpus, .gettext("Cluster")))) * 100,
                 sapply(corpusSubClust$lower, attr, "height"))
    rownames(val) <- c(.gettext("Number of documents"), .gettext("% of documents"), .gettext("Within-cluster variance"))
    colnames(val) <- seq.int(ncol(val))
    names(dimnames(val)) <- c("", .gettext("Cluster"))
    tkinsert(txt, "end", paste(capture.output(format(as.data.frame(val), nsmall=1, digits=2, width=6)),
                               collapse="\n"), "fixed")

    meta <- meta(corpus)[!colnames(meta(corpus)) %in% c("MetaID", .gettext("Cluster"))]
    clusters <- meta(corpus, .gettext("Cluster"))[[1]]

    # Set by createClassesDlg()
    # It is more correct to use exactly the same matrix, and it is more efficient
    if(length(attr(corpusSubClust, "sparsity")) > 0)
        sparsity <- attr(corpusSubClust, "sparsity")
    else
        sparsity <- 1

    clusterDtm <- suppressWarnings(rollup(dtm, 1, clusters))

    if(nterms > 0) {
        specTerms <- specificTerms(dtm, clusters, p, nterms, sparsity, min.occ)
    }

    for(j in 1:ncol(val)) {
        if(nterms > 0) {
            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Terms specific of cluster %i:"), j), "\n", sep=""),
                     "heading")
            tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
            tkinsert(listbox, "end", sprintf(.gettext("Cluster %i:"), j))
            tkitemconfigure(listbox, mark, background="grey")
            mark <- mark + 1

            tkinsert(txt, "end", paste(capture.output(print(specTerms[[j]], nsmall=2, digits=2,
                                                             width=max(nchar(colnames(specTerms[[j]]), "width")))),
                                       collapse="\n"), "fixed")
        }

        if(ndocs > 0) {
            # Remove terms that do not appear in the cluster
            keep <- as.matrix(clusterDtm[j,] > 0)
            docDtm <- dtm[clusters %in% j, keep]
            clustDtm <- clusterDtm[j, keep]
            dev <- sweep(as.matrix(docDtm)/row_sums(docDtm), 2,
                         as.matrix(clustDtm)/sum(clustDtm), "-")
            chisq <- rowSums(sweep(dev^2, 2, col_sums(dtm[,keep])/sum(dtm), "/"))
            chisq <- sort(chisq)[seq.int(1, min(length(chisq), ndocs))]
            docs <- names(chisq)

            tkinsert(txt, "end",
                     paste("\n\n", sprintf(.gettext("Documents specific of cluster %i:"), j), "\n", sep=""),
                     "heading")

            df <- data.frame(row.names=docs, chisq)
            colnames(df) <- .gettext("Chi2 dist. to centroid")

            tkinsert(txt, "end", paste(capture.output(format(df, nsmall=1, digits=2)), collapse="\n"), "fixed")

            tkinsert(txt, "end", "\n", "heading")

            # We need to use IDs rather than indexes to access documents in the corpus
            # since some documents may have been skipped in the clustering
            for(id in docs) {
                tkinsert(txt, "end", paste("\n", id, "\n", sep=""),
                         "articlehead")
                tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
                mark <- mark + 1
                tkinsert(listbox, "end", id)

                origin <- meta(corpus[[id]], "Origin")
                date <- meta(corpus[[id]], "DateTimeStamp")
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
    }

    if(ncol(meta) > 0) {
        tkinsert(txt, "end",
                  paste("\n\n", .gettext("Distribution of variables among clusters:"), "\n", sep=""),
                 "heading")
        tkinsert(txt, "end", paste(.gettext("Row %"), "\n", sep=""), "details")
        tkinsert(listbox, "end", .gettext("Variables"))
        tkitemconfigure(listbox, mark, background="grey")
        tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
        mark <- mark + 1

        # Keep in sync with corpusCa()

        dupLevels <- any(duplicated(unlist(lapply(meta,
            function(x) substr(unique(as.character(x[!is.na(x)])), 0, 30)))))

        # Just in case variables have common levels, and are truncated to the same string
        vars <- colnames(meta)
        vars10<-make.unique(substr(vars, 0, 10))
        vars20<-make.unique(substr(vars, 0, 20))

        tab <- lapply(1:length(vars),
                      function(i) {
                          var <- vars[i]

                          # We call factor() to drop empty levels, if any
                          mat <- with(meta(corpus), table(factor(get(var)), factor(get(.gettext("Cluster")))))

                          # If only one level is present, don't add the level name (e.g. YES),
                          # except if all values are the same (in which case variable is useless but is more obvious that way)
                          if(nrow(mat) == 1 && any(is.na(meta(corpus, .gettext("Cluster"))[[1]])))
                              rownames(mat)<-vars20[i]
                          # In case of ambiguous levels of only numbers in levels, add variable names everywhere
                          else if(dupLevels || !any(is.na(suppressWarnings(as.numeric(rownames(mat))))))
                              rownames(mat)<-make.unique(paste(vars10[i], substr(rownames(mat), 0, 30)))
                          else # Most general case: no need to waste space with variable names
                              rownames(mat)<-substr(rownames(mat), 0, 30)

                          mat
                       })

        tab <- do.call(rbind, tab)
        tab <- rbind(tab, colSums(tab))
        rownames(tab)[nrow(tab)] <- .gettext("Corpus")
        names(dimnames(tab)) <- c("", .gettext("Cluster"))
        tab <- prop.table(tab, 1) * 100

        tkinsert(txt, "end", paste(capture.output(format(as.data.frame(tab), nsmall=1, digits=2, width=5)),
                                   collapse="\n"), "fixed")
    }

    # Only raise the window when we're done, as filling it may take some time
    tkraise(window)

    return()
}

corpusClustDlg <- function() {
    initializeDialog(title=.gettext("Run Hierarchical Clustering"))

    haveCa <- exists("corpusCa")

    setState <- function(...) {
        if(tclvalue(tclType) == "full") {
            caState <- "disabled"
            fullState <- "normal"
        }
        else {
            caState <- "normal"
            fullState <- "disabled"
        }

        tkconfigure(spinSparsity, state=fullState)
        tkconfigure(sparsityLabel, state=fullState)
        tkconfigure(labelNDocs, state=fullState)

        tkconfigure(sliderDim, state=caState)
        tkconfigure(dimLabel, state=caState)
    }

    tclType <- tclVar("full")
    fullButton <- ttkradiobutton(top, variable=tclType,
                                 value="full", text=.gettext("Use full document-term matrix"),
                                 command=setState)
    caButton <- ttkradiobutton(top, variable=tclType,
                               value="ca", text=.gettext("Use dimensions of correspondence analysis"),
                               state=if(haveCa) "active" else "disabled", command=setState)

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

    tclSparsity <- tclVar(100 - ceiling(1/nrow(dtm) * 100))
    spinSparsity <- tkwidget(top, type="spinbox", from=0, to=100,
                             inc=0.1, textvariable=tclSparsity,
                             validate="all", validatecommand=function(P) .validate.unum(P, fun=updateNDocs))

    sparsityLabel <- labelRcmdr(top, text=.gettext("Remove terms missing from more than (% of documents):"))
    updateNDocs(tclvalue(tclSparsity))


    tclDim <- tclVar(2)
    sliderDim <- tkscale(top, from=1, to=if(haveCa) corpusCa$nd else 5,
                         showvalue=TRUE, variable=tclDim,
		         resolution=1, orient="horizontal",
                         state="disabled")
    dimLabel <- labelRcmdr(top, text=.gettext("Retain dimensions from 1 to:"), state="disabled")


    onOK <- function() {
        type <- tclvalue(tclType)
        sparsity <- as.numeric(tclvalue(tclSparsity))/100
        dim <- tclvalue(tclDim)

        if(is.na(sparsity) || sparsity <= 0 || sparsity > 1) {
            .Message(.gettext("Please specify a sparsity value between 0 (excluded) and 100%."), type="error")
            return()
        }

        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        # removeSparseTerms() does not accept 1
        if(type == "ca") {
            if(is.null(corpusCa$rowvars))
                doItAndPrint(sprintf('chisqDist <- dist(sweep(corpusCa$rowcoord[, 1:%s, drop=FALSE], 2, corpusCa$sv[1:%s], "*"))',
                                     dim, dim))
            else
                doItAndPrint(sprintf('chisqDist <- dist(sweep(corpusCa$rowcoord[-corpusCa$rowvars, 1:%s, drop=FALSE], 2, corpusCa$sv[1:%s], "*"))',
                                     dim, dim))

            # Memory allocation can fail here, so avoid stacking errors
            if(!exists("chisqDist"))
                return()

            doItAndPrint('corpusClust <- hclust(chisqDist, method="ward")')
            doItAndPrint(sprintf('attr(corpusClust, "caDim") <- %s', dim))

            doItAndPrint(sprintf('attr(corpusClust, "nterms") <- %s', nrow(corpusCa$colcoord)))

            doItAndPrint("rm(chisqDist)")
            gc()
        }
        else if(sparsity < 1) {
            doItAndPrint(sprintf("clustDtm <- removeSparseTerms(dtm, %s)", sparsity))

            if(any(row_sums(clustDtm) == 0)) {
                Message(paste(.gettext("Documents skipped from hierarchical clustering:\n"),
                              paste(rownames(clustDtm)[row_sums(clustDtm) == 0], collapse=", ")),
                        type="note")

		    .Message(sprintf(.gettext("%i documents have been skipped because they do not include any occurrence of the terms retained in the final document-term matrix. Their list is available in the \"Messages\" area.\n\nIncrease the value of the 'sparsity' parameter if you want to include them."),
		                     sum(row_sums(clustDtm) == 0)),
		             type="info")

                doItAndPrint('clustDtm <- clustDtm[row_sums(clustDtm) > 0,]')
            }

            doItAndPrint('chisqDist <- dist(sweep(clustDtm/row_sums(clustDtm), 2, sqrt(col_sums(clustDtm)/sum(clustDtm)), "/"))')

            # Memory allocation can fail here, so avoid stacking errors
            if(!exists("chisqDist")) {
                doItAndPrint("rm(clustDtm)")
                return()
            }

            doItAndPrint('corpusClust <- hclust(chisqDist, method="ward")')

            # Used by createClustersDialog() and showCorpusClust() to recreate the dtm
            doItAndPrint(sprintf('attr(corpusClust, "sparsity") <- %s', sparsity))

            doItAndPrint(sprintf('attr(corpusClust, "nterms") <- %s', ncol(clustDtm)))

            doItAndPrint("rm(clustDtm, chisqDist)")

            gc()
        }
        else {
            doItAndPrint('chisqDist <- dist(sweep(dtm/row_sums(dtm), 2, sqrt(col_sums(dtm)/sum(dtm)), "/"))')

            # Memory allocation can fail here, so avoid stacking errors
            if(!exists("chisqDist"))
                return()

            doItAndPrint('corpusClust <- hclust(chisqDist, method="ward")')
            doItAndPrint(sprintf('attr(corpusClust, "nterms") <- %s', ncol(dtm)))
            doItAndPrint("rm(chisqDist)")
            gc()
        }

        # For the Create clusters item
        activateMenus()

        tkfocus(CommanderWindow())

        createClustersDlg()
    }

    OKCancelHelp(helpSubject="corpusClustDlg")
    tkgrid(fullButton, sticky="sw")
    tkgrid(sparsityLabel, spinSparsity, sticky="sw", pady=c(0, 6), padx=c(24, 0))
    tkgrid(labelNDocs, sticky="sw", pady=6, columnspan=2, padx=c(24, 0))
    tkgrid(caButton, sticky="sw", pady=c(12, 0))
    tkgrid(dimLabel, sliderDim, sticky="sw", pady=c(0, 6), padx=c(24, 0))
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix()
}

createClustersDlg <- function(..., plot=TRUE) {
    if(!(exists("corpusClust") && class(corpusClust) == "hclust")) {
        .Message(message=.gettext("Please run a hierarchical clustering on the corpus first."),
                type="error")
        return()
    }

    if(plot) {
        setBusyCursor()

        # Do not plot all leafs if there are too many of them (can even crash!)
        if(length(corpusClust$labels) > 500) {
            # If 500th is already at 0, cutting there won't have any effect,
            # so take the last non-0 height in the 1:500 range
            revList <- rev(corpusClust$height)[1:500]
            height <- floor(tail(revList[revList > sqrt(.Machine$double.eps)], 1) * 1e4)/1e4
            doItAndPrint(sprintf('plot(cut(as.dendrogram(corpusClust), h=%s)$upper, leaflab="none", ylab="%s", main="%s")',
                                 height,
                                 .gettext("Within-cluster variance"),
                                 .gettext("Upper part of documents dendrogram")))
        }
        else {
            doItAndPrint(sprintf('plot(as.dendrogram(corpusClust), nodePar=list(pch=NA, lab.cex=0.8), %sylab="%s", main="%s")',
                                 if(length(corpus) > 20) 'leaflab="none", ' else "",
                                 .gettext("Within-cluster variance"),
                                 .gettext("Full documents dendrogram")))
        }

        setIdleCursor()
    }

    initializeDialog(title=.gettext("Create Clusters"))

    tclNClust <- tclVar(5)
    sliderNClust <- tkscale(top, from=2, to=min(15, length(corpusClust$order)),
                            showvalue=TRUE, variable=tclNClust,
		                    resolution=1, orient="horizontal")

    tclNDocs <- tclVar(5)
    spinNDocs <- tkwidget(top, type="spinbox", from=1, to=length(corpus),
                          inc=1, textvariable=tclNDocs,
                          validate="all", validatecommand=.validate.uint)

    tclNTerms <- tclVar(20)
    spinNTerms <- tkwidget(top, type="spinbox", from=1, to=ncol(dtm),
                           inc=1, textvariable=tclNTerms,
                           validate="all", validatecommand=.validate.uint)

    tclP <- tclVar(10)
    spinP <- tkwidget(top, type="spinbox", from=0, to=100,
                      inc=1, textvariable=tclP,
                      validate="all", validatecommand=.validate.unum)

    tclOcc <- tclVar(2)
    spinOcc <- tkwidget(top, type="spinbox", from=1, to=.Machine$integer.max,
                        inc=1, textvariable=tclOcc,
                        validate="all", validatecommand=.validate.uint)

    onOK <- function() {
        closeDialog()

        setBusyCursor()
        on.exit(setIdleCursor())

        nclust <- as.numeric(tclvalue(tclNClust))
        ndocs <- as.numeric(tclvalue(tclNDocs))
        nterms <- as.numeric(tclvalue(tclNTerms))
        p <- as.numeric(tclvalue(tclP))/100
        occ <- as.numeric(tclvalue(tclOcc))
        height <- floor(rev(corpusClust$height)[nclust-1] * 1e4)/1e4

        doItAndPrint(paste("corpusSubClust <- cut(as.dendrogram(corpusClust), h=",
                           height, ")", sep=""))

        # Hack needed to replace cutree(), which orders clusters by order of appearance in the data,
        # which does not correspond to that used by the dendrogram
        doItAndPrint(sprintf('clusters <- rep(1:%s, sapply(corpusSubClust$lower, attr, "members"))', nclust))
        doItAndPrint('names(clusters) <- rapply(corpusSubClust$lower, function(x) attr(x, "label"))')

        # If some documents were skipped, we need to fill with NA
        doItAndPrint(sprintf('meta(corpus, "%s") <- NA', .gettext("Cluster")))

        doItAndPrint(sprintf('meta(corpus, "%s")[match(names(clusters), names(corpus), nomatch=0),] <- clusters',
                             .gettext("Cluster")))

        # If corpus was split, we cannot add cluster back into corpusVars
        if(exists("corpusVars")) {
            # If corpus was split, we cannot add cluster back into corpusVars
            if(nrow(corpusVars) == length(corpus))
                doItAndPrint(paste("corpusVars$",  .gettext("Cluster"),
                                   " <- meta(corpus, tag=\"", .gettext("Cluster"), "\")[[1]]", sep=""))
        }
        else {
            doItAndPrint(paste("corpusVars <- data.frame(",  .gettext("Cluster"),
                               "=meta(corpus, tag=\"", .gettext("Cluster"), "\")[[1]])", sep=""))
        }

        # Set by corpusClustDlg() and used by showCorpusClustering() to recreate dtm
        if(length(attr(corpusClust, "sparsity")) > 0)
            doItAndPrint(sprintf('attr(corpusSubClust, "sparsity") <- %s', attr(corpusClust, "sparsity")))

        doItAndPrint(sprintf('attr(corpusSubClust, "nterms") <- %s', attr(corpusClust, "nterms")))

        if(length(attr(corpusClust, "caDim")) > 0)
            doItAndPrint(sprintf('attr(corpusSubClust, "caDim") <- %s', attr(corpusClust, "caDim")))

        doItAndPrint(sprintf('plot(corpusSubClust$upper, nodePar=list(pch=NA, lab.cex=0.8), ylab="%s", main="%s")',
                             .gettext("Within-cluster variance"),
                             .gettext("Clusters dendrogram")))
        doItAndPrint(sprintf("showCorpusClustering(corpusSubClust, %i, %i, %f, %i)", ndocs, nterms, p, occ))
        doItAndPrint("rm(clusters)")

        tkfocus(CommanderWindow())
    }

    OKCancelHelp(helpSubject="createClustersDlg")
    tkgrid(.titleLabel(top, text=.gettext("Clusters creation:")),
           sticky="sw", pady=0)
    tkgrid(labelRcmdr(top, text=.gettext("Number of clusters to retain:")), sliderNClust,
           sticky="sw", pady=c(0, 6), padx=c(6, 0))
    tkgrid(.titleLabel(top, text=.gettext("Documents specific of clusters:")),
           sticky="sw", pady=c(24, 0))
    tkgrid(labelRcmdr(top, text=.gettext("Maximum number of documents to show per cluster:")), spinNDocs,
           sticky="sw", pady=c(0, 6), padx=c(6, 0))
    tkgrid(.titleLabel(top, text=.gettext("Terms specific of clusters:")),
           sticky="sw", pady=c(24, 0))
    tkgrid(labelRcmdr(top, text=.gettext("Show terms with a probability below (%):")), spinP,
           sticky="sw", pady=6, padx=c(6, 0))
    tkgrid(labelRcmdr(top, text=.gettext("Only retain terms with a number of occurrences above:")), spinOcc,
           sticky="sw", pady=6, padx=c(6, 0))
    tkgrid(labelRcmdr(top, text=.gettext("Maximum number of terms to show per cluster:")), spinNTerms,
           sticky="sw", pady=6, padx=c(6, 0))
    tkgrid(buttonsFrame, columnspan=2, sticky="ew", pady=6)
    dialogSuffix()
}

