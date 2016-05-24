`truthTable` <-
function(data, outcome = "", conditions = "", n.cut = 1,
         incl.cut = 1, complete = FALSE, show.cases = FALSE,
         sort.by = "", use.letters = FALSE, inf.test = "", ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    metacall <- match.call(expand.dots = TRUE)
    other.args <- list(...)
    via.pof <- "via.pof" %in% names(other.args)
    
    ica <- 1
    
    if (is.character(incl.cut) & length(incl.cut) == 1) {
        incl.cut <- QCA::splitstr(incl.cut)
    }
    
    icp <- incl.cut[1]
    if (length(incl.cut) > 1) {
        ica <- incl.cut[2]
    }
    
    
    ### 
    ### ### backwards compatibility 
    ### 
        neg.out <- FALSE
        if ("neg.out" %in% names(other.args)) {
            neg.out <- other.args$neg.out
        }
        
        if ("incl.cut1" %in% names(other.args) & identical(icp, 1)) {
            icp <- other.args$incl.cut1
            incl.cut[1] <- icp
        }
        
        if ("incl.cut0" %in% names(other.args) & identical(ica, 1)) {
            ica <- other.args$incl.cut0
            incl.cut[2] <- ica
        }
    ### 
    ### ### backwards compatibility 
    ### 
    
    names(data) <- toupper(names(data))
    conditions <- toupper(conditions)
    outcome <- toupper(outcome)
    
    if (length(outcome) > 1) {
        cat("\n")
        stop(simpleError("Only one outcome is allowed.\n\n"))
    }
    
    outcome.copy <- outcome
    initial.data <- data
    
    if (substring(outcome, 1, 1) == "~") {
        neg.out <- TRUE
        outcome <- substring(outcome, 2)
    }
    
    # for the moment, toupper(outcome) is redundant but should the negation be 
    # treated with lower case letters in the future, it will prove important
    if (!identical(outcome, "")) {
        if (! toupper(QCA::curlyBrackets(outcome, outside=TRUE)) %in% colnames(data)) {
            cat("\n")
            stop(simpleError("Inexisting outcome name.\n\n"))
        }
    }
    
    if (grepl("\\{|\\}", outcome)) {
        outcome.value <- QCA::curlyBrackets(outcome)
        outcome <- QCA::curlyBrackets(outcome, outside=TRUE)
        
        data[, toupper(outcome)] <- as.numeric(data[, toupper(outcome)] %in% QCA::splitstr(outcome.value))
    }
    ### this was supposed to treat the negation using lower case letters
    # else if (! outcome %in% colnames(data)) {
    #     data[, toupper(outcome)] <- 1 - data[, toupper(outcome)]
    # }
    
    # already on line 23
    # outcome <- toupper(outcome)
    ### 
    
    if (identical(conditions, "")) {
        conditions <- names(data)[-which(names(data) == outcome)]
    }
    else {
        if (is.character(conditions) & length(conditions) == 1) {
            conditions <- QCA::splitstr(conditions)
        }
    }
    
    if (is.character(sort.by) & length(sort.by) == 1 & !identical(sort.by, "")) {
        sort.by <- QCA::splitstr(sort.by)
    }
    
    decreasing <- TRUE # just to set a default value
    if ("decreasing" %in% names(other.args)) {
        decreasing <- other.args$decreasing
    }
    
    if (is.character(decreasing) & length(decreasing) == 1) {
        decreasing <- QCA::splitstr(decreasing)
    }
    
    
    if (!identical(inf.test, "")) {
        inf.test <- QCA::splitstr(inf.test)
    }
    
    if (!via.pof) {
        verify.tt(data, outcome, conditions, complete, show.cases, icp, ica, inf.test)
    }
    
    data <- data[, c(conditions, outcome)]
    
    if (ica > icp) {
        ica <- icp
    }
    
    colnames(data) <- toupper(colnames(data))
    colnames(initial.data) <- toupper(colnames(initial.data))
    conditions <- toupper(conditions)
    outcome <- toupper(outcome)
    
    initial.data <- initial.data[, c(conditions, outcome)]
    
    
    if (neg.out) {
        data[, outcome] <- 1 - data[, outcome]
    }
    
    dc.code <- unique(unlist(lapply(data, function(x) {
        if (is.numeric(x)) {
            return(x[x < 0])
        }
        else {
            return(as.character(x[x %in% c("-", "dc")]))
        }
    })))
    
    if (length(dc.code) == 0) {
        dc.code <- -1
    }
    else if (length(dc.code) > 1) {
        cat("\n")
        stop(simpleError("Multiple \"Don't care\" codes found.\n\n"))
    }
    
    data <- as.data.frame(lapply(data, function(x) {
        # to make sure that any factor is character
        x <- as.character(x)
        
        # otherwise replacement was not possible
        x[x == dc.code] <- -1
        
        return(QCA::asNumeric(x))
    }))
    
    names(data) <- c(conditions, outcome)
    
    data[data < 0] <- -1
    rownames(data) <- rownames(initial.data)
    
    nofconditions <- length(conditions)
    # fuzzy.cc <- apply(data[, conditions, drop=FALSE], 2, function(x) any(abs(x - round(x)) >= .Machine$double.eps^0.5))
    
    fuzzy.cc <- apply(data[, conditions, drop=FALSE], 2, function(x) any(x %% 1 > 0))
    
    
    for (i in seq(length(conditions))) {
        if (!fuzzy.cc[i]) {
            copy.cc <- data[, i]
            if (any(copy.cc < 0)) {
                copy.cc[copy.cc < 0] <- max(copy.cc) + 1
                data[, i] <- copy.cc
            }
        }
    }
    
    # the data MUST begin with 0 and MUST be incremented by 1 for each level...!
    # perhaps trying something like
    # apply(data[, conditions], 2, function(x) length(unique(x))) + 1
    noflevels <- apply(data[, conditions, drop=FALSE], 2, max) + 1
    noflevels[noflevels == 1] <- 2
    noflevels[fuzzy.cc] <- 2
    
    
    if (via.pof) {
        return(as.vector(noflevels))
    }
    
    
    condata <- data[, conditions]
    if (any(fuzzy.cc)) {
        condata[, fuzzy.cc] <- lapply(condata[, fuzzy.cc], function(x) {
            (x > 0.5)*1
        })
    }
    
    line.data <- as.vector(as.matrix(condata) %*% c(rev(cumprod(rev(noflevels))), 1)[-1])
    condata <- condata[order(line.data), ]
    uniq <- which(!duplicated(condata))
    tt <- condata[uniq, ]
    
    rownstt <- sort(line.data)[uniq] + 1
    rownames(tt) <- rownstt
    
    # ipc <- .Call("truthTable", as.matrix(data[, conditions]), as.matrix(tt), as.numeric(fuzzy.cc), data[, outcome], PACKAGE="QCAGUI")
    ipc <- QCA::callTruthTable(as.matrix(data[, conditions]), as.matrix(tt), as.numeric(fuzzy.cc), data[, outcome])
    colnames(ipc) <- rownstt
    
    exclude <- ipc[1, ] < n.cut
    
    if (sum(!exclude) == 0) {
        cat("\n")
        stop(simpleError("There are no combinations at this frequency cutoff.\n\n"))
    }
    
    tt$OUT <- "?"
    tt$OUT[!exclude] <- as.numeric(ipc[2, ] >= (icp - .Machine$double.eps ^ 0.5))[!exclude]
    tt$OUT[!exclude & ipc[2, ] < icp & ipc[2, ] >= (ica - .Machine$double.eps ^ 0.5)] <- "C"
    tt <- cbind(tt, t(ipc))
    
    cases <- sapply(line.data[order(line.data)][uniq], function(x) {
        paste(rownames(data)[which(line.data == x)], collapse=",")
    })
    
    excluded <- tt[exclude, , drop = FALSE]
    
    if (length(conditions) < 8) {
        
        ttc <- as.data.frame(matrix(nrow = prod(noflevels), ncol = ncol(tt)))
        colnames(ttc) <- colnames(tt)
        ttc[, seq(length(conditions))] <- createMatrix(noflevels)
        ttc$OUT   <- "?"
        ttc$n     <-  0
        ttc$incl  <- "-"
        
        # sometimes a causal condition may be named PRI (see Porter data)
        whichpri <- which(colnames(ttc) == "PRI")
        
        ttc[, whichpri[length(whichpri)]] <- "-"  
        ttc[rownames(tt), ] <- tt
        tt <- ttc
    }
    
    
    if (!identical(sort.by, "")) {
        
        if (is.logical(sort.by)) { # & !is.null(names(sort.by)) # if logical, it should _always_ have names
            decreasing <- as.vector(sort.by)
            sort.by <- names(sort.by)
        }
        else {
            # just to make sure we _do_ have a "decreasing" object
            if (missing(decreasing)) {
                decreasing <- rep(TRUE, length(sort.by))
            }
            else {
                if (is.logical(decreasing)) {
                    if (length(decreasing) == 1) {
                        decreasing <- rep(decreasing, length(sort.by))
                    }
                    else if (length(decreasing) < length(sort.by)) {
                        decreasing <- c(decreasing, rep(TRUE, length(sort.by) - length(decreasing)))
                    }
                }
                else {
                    decreasing <- rep(TRUE, length(sort.by))
                }
            }
        }
        
        sort.by[sort.by == "out"] <- "OUT"
        
        decreasing <- decreasing[sort.by %in% names(tt)]
        sort.by <- sort.by[sort.by %in% names(tt)]
        
        
        rowsorder <- seq_len(nrow(tt))
        for (i in rev(seq(length(sort.by)))) {
            rowsorder <- rowsorder[order(tt[rowsorder, sort.by[i]], decreasing = decreasing[i])]
        }
        
        sortvector <- rep(1, nrow(tt))
        sortvector[tt[rowsorder, "OUT"] == "?"] <- 2
        rowsorder <- rowsorder[order(sortvector)]
        
    }
    
    
    
    for (i in seq(length(conditions))) {
        if (!fuzzy.cc[i]) {
            if (any(initial.data[, i] == dc.code)) {
                tt[, i][tt[, i] == max(tt[, i])] <- dc.code
                data[, i][data[, i] == max(data[, i])] <- dc.code
                noflevels[i] <- noflevels[i] - 1
            }
        }
    }
    
    statistical.testing <- FALSE
    
    if (inf.test[1] == "binom") {
        statistical.testing <- TRUE
        if (length(inf.test) > 1) {
            alpha <- as.numeric(inf.test[2]) # already checked if a number between 0 and 1
        }
        else {
            alpha <- 0.05
        }
        
        observed <- which(tt$OUT != "?")
        success <- round(tt[observed, "n"] * as.numeric(tt[observed, "incl"]))
        
        tt$pval1 <- "-"
        if (length(incl.cut) > 1) {
            tt$pval0 <- "-"
        }
        tt[observed, "OUT"] <- 0
        
        for (i in seq(length(observed))) {
            
            pval1 <- tt[observed[i], "pval1"] <- binom.test(success[i], tt[observed[i], "n"], p = icp, alternative = "greater")$p.value
            if (length(incl.cut) > 1) {
                pval0 <- tt[observed[i], "pval0"] <- binom.test(success[i], tt[observed[i], "n"], p = ica, alternative = "greater")$p.value
            }
            
            if (pval1 < alpha) {
                tt[observed[i], "OUT"] <- 1
            }
            else if (length(incl.cut) > 1) {
                if (pval0 < alpha) {
                    tt[observed[i], "OUT"] <- "C"
                }
            }
        }
    }
    
    
    # deal with the show.cases in the print function
    # if (show.cases) {
    
        # apparently this is necessary, otherwise the new column cases will be a factor
        tt$cases <- ""
        
        if (length(conditions) < 8) {
            tt$cases[rownstt] <- cases
        }
        else {
            tt$cases <- cases[!exclude]
        }
        
    # }
    
    
    x <- list(tt = tt, indexes = rownstt, noflevels = as.vector(noflevels),
              initial.data = initial.data, recoded.data = data, cases = cases, 
              options = list(outcome = outcome.copy, conditions = conditions, neg.out = neg.out, n.cut = n.cut,
                             incl.cut = incl.cut, complete = complete, show.cases = show.cases,
                             use.letters = use.letters, inf.test = statistical.testing,
                             incl.cut1 = incl.cut[1], incl.cut0 = ifelse(length(incl.cut) == 2, incl.cut[2], incl.cut)))
    
    if (any(exclude)) {
        excluded$cases <- ""
        excluded$cases <- cases[exclude]
        x$excluded <- structure(list(tt = excluded,
                                    options = list(show.cases = TRUE, complete = FALSE, excluded = TRUE)), class="tt")
    }
    
    # also verify if not already letters
    if (use.letters & any(nchar(conditions) > 1)) { 
        colnames(x$tt)[seq(nofconditions)] <- LETTERS[seq(nofconditions)]
    }
    
    # PRI <- FALSE
    # if ("PRI" %in% names(other.args)) {
    #     if (is.logical(other.args$PRI)) {
    #         PRI <- other.args$PRI[1]
    #     }
    # }
    
    # x$options$PRI <- PRI
    
    if (!identical(sort.by, "")) {
        x$rowsorder <- rowsorder
    }
    
    x$origin <- "QCAGUI"
    x$call <- metacall
    return(structure(x, class="tt"))
}

