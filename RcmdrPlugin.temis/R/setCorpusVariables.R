doSetCorpusVariables <- function() doItAndPrint("setCorpusVariables()")

setCorpusVariables <- function() {
    if(!exists("corpus") || !("Corpus" %in% class(corpus))) {
        .Message(message=.gettext("Please import a corpus first."),
                 type="error")
        return()
    }

    if(!activeDataSetP()) {
        .Message(message=.gettext("Please create or import a data set first."),
                 type="error")
        return()
    }

    if(!checkVariables()) {
        .Message(message=.gettext("Please create at least one variable (column)."),
                 type="error")
        return()
    }

    # If corpus was split, we need to replicate variables
    split <- isTRUE(meta(corpus, type="corpus", tag="split"))

    dset <- get(ActiveDataSet())
    len <- if (split) length(unique(meta(corpus, .gettext("Doc N"))[[1]])) else length(corpus)
    if(nrow(dset) != len) {
        .Message(message=sprintf(.gettext("Active data set must contain exactly %d rows."), len),
                 type="error")
        return()
    }

    # Remove dropped and empty variables
    for(var in colnames(meta(corpus))[!colnames(meta(corpus)) %in%
            c(colnames(dset), .gettext("Doc N"), .gettext("Doc ID"), .gettext("Cluster"),
              sapply(corpusVars, function(x) all(is.na(x) | x == "")))])
        meta(corpus, var) <<- NULL

    # Add new variables
    indices <- which(sapply(dset, function(x) !all(is.na(x) | x == "", na.rm=TRUE)))

    # Drop empty levels, which notably happen when changing values manually
    if(length(indices) > 0) {
        if(split) {
            for(i in indices)
               meta(corpus, colnames(dset)[i]) <<- droplevels(factor(dset[meta(corpus, .gettext("Doc N"))[[1]], i]))
        }
        else {
            for(i in indices)
                meta(corpus, colnames(dset)[i]) <<- droplevels(factor(dset[[i]]))
        }
    }

    # Update names only if they changed
    oldDocNames <- if(split) unique(meta(corpus, .gettext("Doc ID"))[[1]]) else names(corpus)
    corpusNames <- names(corpus)
    if(!identical(oldDocNames, row.names(dset))) {
        if(split) {
            names(corpus) <<- make.unique(row.names(dset)[meta(corpus, .gettext("Doc N"))[[1]]])
            meta(corpus, .gettext("Doc ID")) <<- row.names(dset)[meta(corpus, .gettext("Doc N"))[[1]]]
        }
        else {
            names(corpus) <- row.names(dset)
        }

        # Update the names of the dtm since it affects all operations and cannot be done manually
        # We assume the dtm corresponds to the current corpus if names were identical
        if(identical(corpusNames, rownames(dtm)))
            justDoIt("rownames(dtm) <- names(corpus)")

        if(exists("wordsDtm") && identical(corpusNames, rownames(wordsDtm)))
            justDoIt("rownames(wordsDtm) <- names(corpus)")
    }
}
