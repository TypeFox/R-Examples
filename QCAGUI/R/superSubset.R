`superSubset` <-
function(data, outcome = "", conditions = "", relation = "nec", incl.cut = 1,
    cov.cut = 0, use.tilde = FALSE, use.letters = FALSE, ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    memcare <- FALSE # to be updated with a future version
    
    other.args <- list(...)
    
    colnames(data) <- toupper(colnames(data))
    
    ### backwards compatibility 
        neg.out <- FALSE
        if ("neg.out" %in% names(other.args)) {
            neg.out <- other.args$neg.out
        }
    ### 
    
    
    
    incl.cut <- incl.cut - .Machine$double.eps ^ 0.5
    if (cov.cut > 0) {
        cov.cut <- cov.cut - .Machine$double.eps ^ 0.5
    }
    
    ### TO MODIFY when the negation will be treated with lower case letters
    outcome <- toupper(outcome)
    
    if (substring(outcome, 1, 1) == "~") {
        neg.out <- TRUE
        outcome <- substring(outcome, 2)
    }
    
    # for the moment, toupper(outcome) is redundant but if in the future
    # the negation will be treated with lower case letters, it will prove important
    if (! toupper(QCA::curlyBrackets(outcome, outside=TRUE)) %in% toupper(colnames(data))) {
        cat("\n")
        stop(simpleError("Inexisting outcome name.\n\n"))
    }
    
    if (grepl("\\{|\\}", outcome)) {
        outcome.value <- QCA::curlyBrackets(outcome)
        outcome <- QCA::curlyBrackets(outcome, outside=TRUE)
        
        data[, toupper(outcome)] <- as.numeric(data[, toupper(outcome)] %in% QCA::splitstr(outcome.value))
    }
    ### this is supposed to treat the negation using lower case letters, maybe in a future version
    # else if (! outcome %in% colnames(data)) {
    #     data[, toupper(outcome)] <- 1 - data[, toupper(outcome)]
    # }
    
    # to undo line 45
    # outcome <- toupper(outcome)
    ### 
    
    if (identical(conditions, "")) {
        conditions <- names(data)[-which(names(data) == outcome)]
    }
    else {
        conditions <- QCA::splitstr(conditions)
    }
    
    conditions <- toupper(conditions)
    
    verify.data(data, outcome, conditions)
    
    if (!(QCA::nec(relation) | QCA::suf(relation) | relation %in% c("sufnec", "necsuf"))) {
        stop(simpleError("\nThe relationship should be \"necessity\", \"sufficiency\", \"sufnec\" or \"necsuf\".\n\n"))
    }
    
    relationcopy <- relation
    if (relation %in% c("sufnec", "necsuf")) {
        cov.cut <- incl.cut
    }
    
    if (relation == "sufnec") {
        relation <- "sufficiency"
    }
    else if (relation == "necsuf") {
        relation <- "necessity"
    }
    
    colnames(data) <- toupper(colnames(data))
    conditions <- replacements <- toupper(conditions)
    outcome <- toupper(outcome)
    
    data <- data[, c(conditions, outcome)]
    nofconditions <- length(conditions)
    
    
    if (neg.out) {
        data[, outcome] <- 1 - data[, outcome]
    }
    
    uplow <- !use.tilde
    
    fc <- apply(data[, conditions], 2, function(x) any(x %% 1 > 0))
    
    if (mv.data <- any(data[, conditions] > 1)) {
        uplow <- use.tilde <- FALSE
    }
    
    alreadyletters <- sum(nchar(conditions)) == length(conditions)
    
    collapse <- ifelse(alreadyletters & uplow | use.tilde, "", "*")
    
    if (use.letters & !alreadyletters) {
        replacements <- LETTERS[seq(length(conditions))]
        names(replacements) <- conditions
        colnames(data)[seq(length(conditions))] <- conditions <- replacements
        collapse <- ifelse(!uplow | use.tilde, "*", "")
    }
    
    noflevels <- apply(data[, conditions], 2, max) + 1L
    noflevels[fc] <- 2
    mbase <- c(rev(cumprod(rev(noflevels + 1L))), 1)[-1]
    
    
    if (memcare) {
        # CMatrix <- .Call("superSubsetMem", as.matrix(data[, conditions]), noflevels, mbase, as.numeric(fc), data[, outcome], as.numeric(relation == "necessity"), PACKAGE="QCAGUI")
        CMatrix <- QCA::callSuperSubsetMem(as.matrix(data[, conditions]), noflevels, mbase, as.numeric(fc), data[, outcome], as.numeric(QCA::nec(relation)))
    }
    else {
        nk <- createMatrix(noflevels + 1L)
        # colnames(nk) <- conditions
        nk <- nk[-1, ] # first row is always empty
        
        # CMatrix <- .Call("superSubset", as.matrix(data[, conditions]), nk, as.numeric(fc), data[, outcome], as.numeric(relation == "necessity"), PACKAGE="QCAGUI")
        CMatrix <- QCA::callSuperSubset(as.matrix(data[, conditions]), nk, as.numeric(fc), data[, outcome], as.numeric(QCA::nec(relation)))
    }
    
    
    # to modify this, attributing colnames copies the object and uses too much memory
    colnames(CMatrix) <- expressions <- seq_len(ncol(CMatrix)) + 1L
    lincl <- ifelse(QCA::nec(relation), 2, 1)
    
    expressions <- expressions[CMatrix[lincl, ] >= incl.cut & CMatrix[3 - lincl, ] >= cov.cut]
    
    prev.result <- FALSE
    lexpressions <- length(expressions)
    
    if (lexpressions > 0) {
        if (QCA::suf(relation)) {
            # expressions <- .Call("removeRedundants", expressions, noflevels, mbase, PACKAGE="QCAGUI")
            expressions <- QCA::callRemoveRedundants(expressions, noflevels, mbase)
        }
        
        result.matrix <- getRow(noflevels + 1L, expressions)
        rownames(result.matrix) <- expressions
        colnames(result.matrix) <- conditions
        result.matrix <- QCA::sortMatrix(result.matrix)
        sum.zeros <- apply(result.matrix, 1, function(idx) sum(idx == 0))
        result.matrix <- result.matrix[order(sum.zeros, decreasing=TRUE), , drop=FALSE]
        row_names <- QCA::writePrimeimp(result.matrix, collapse=collapse, uplow=uplow, use.tilde=use.tilde)
        prev.result <- TRUE
        result <- data.frame(incl  = CMatrix[lincl, rownames(result.matrix)],
             PRI   = CMatrix[5, rownames(result.matrix)],
             cov.r = CMatrix[3 - lincl, rownames(result.matrix)],
             stringsAsFactors=FALSE,
             row.names=row_names)
    }
    
    lexprnec <- 0
    if (QCA::nec(relation)) {
        exprnec <- seq_len(ncol(CMatrix)) + 1L
        
        exprnec <- exprnec[CMatrix[4, ] >= incl.cut & CMatrix[3, ] >= cov.cut]
        
        # exprnec <- .Call("removeRedundants", exprnec, noflevels, mbase, PACKAGE="QCAGUI")
        exprnec <- QCA::callRemoveRedundants(exprnec, noflevels, mbase)
        
        exprnec <- setdiff(exprnec, expressions)
        lexprnec <- length(exprnec)
        
        if (lexprnec + lexpressions == 0) {
            cat("\n")
            stop(simpleError(paste("\nThere are no combinations with incl.cut = ", round(incl.cut, 3), " and cov.cut = ", round(cov.cut, 3), "\n\n", sep="")))
        }
        
        if (lexprnec > 0) {
            result.matrix2 <- getRow(noflevels + 1, exprnec)
            rownames(result.matrix2) <- exprnec
            colnames(result.matrix2) <- conditions
            result.matrix2 <- QCA::sortMatrix(result.matrix2)
            
            sum.zeros <- apply(result.matrix2, 1, function(idx) sum(idx == 0))
            result.matrix2 <- result.matrix2[order(sum.zeros, decreasing=TRUE), , drop=FALSE]
            row_names2 <- QCA::writePrimeimp(result.matrix2, collapse="+", uplow=uplow, use.tilde=use.tilde)
            
            if (prev.result) {
                result <- rbind(result, data.frame(incl  = CMatrix[4, rownames(result.matrix2)],
                    PRI   = CMatrix[6, rownames(result.matrix2)],
                    cov.r = CMatrix[3, rownames(result.matrix2)],
                    stringsAsFactors=FALSE,
                    row.names=row_names2))
                row_names <- c(row_names, row_names2)
                result.matrix <- rbind(result.matrix, result.matrix2)
            }
            else {
                result <- data.frame(incl = CMatrix[4, rownames(result.matrix2)],
                    PRI = CMatrix[6, rownames(result.matrix2)],
                    cov.r = CMatrix[3, rownames(result.matrix2)],
                    stringsAsFactors=FALSE,
                    row.names=row_names2)
                row_names <- row_names2
                result.matrix <- result.matrix2
            }
            
        }
    }
    
    if (lexprnec + lexpressions == 0) {
        cat("\n")
        stop(simpleError(paste("\nThere are no combinations with incl.cut = ", round(incl.cut, 3), " and cov.cut = ", round(cov.cut, 3), "\n\n", sep="")))
    }
    
    
    mins <- matrix(NA, nrow=nrow(data), ncol=nrow(result.matrix))
    for (i in seq(nrow(result.matrix))) {
        mins[, i] <- apply(data[, conditions], 1, function(v) {
            e <- result.matrix[i, , drop=FALSE]
            if (any(ox <- e[fc] == 1)) {
                v[fc][ox] <- 1 - v[fc][ox]
            }
            
            if (length(cp <- v[!fc]) > 0) {
                v[!fc][e[!fc] != cp + 1] <- 0
                v[!fc][e[!fc] == cp + 1] <- 1
            }
            if (rownames(e) %in% expressions) {
                return(min(v[e != 0]))
            }
            else {
                return(max(v[e != 0]))
            }
        })
    }
    
    colnames(mins) <- rownames(result)
    rownames(mins) <- rownames(data)
    mins <- as.data.frame(mins)
    
    if (relationcopy == "sufnec") {
        colnames(result) <- c("inclS", "PRI", "inclN")
    }
    else if (relationcopy == "necsuf") {
        colnames(result) <- c("inclN", "PRI", "inclS")
    }
    
    
    if (QCA::nec(relation)) {
        colnames(result)[2] <- "RoN"
        result[, 2] <- pof(mins, data[, outcome])$incl.cov[, 2]
    }
    
    out.list <- list(incl.cov=result, coms=mins, use.letters=use.letters)
    
    if (use.letters & !alreadyletters) {
        out.list$letters=replacements
    }
    
    out.list$options <- list(
        outcome = outcome,
        neg.out = neg.out,
        conditions = conditions,
        relation = relation,
        incl.cut = incl.cut,
        cov.cut = cov.cut,
        use.tilde = use.tilde,
        use.letters = use.letters
    )
    
    return(structure(out.list, class="sS"))
}



