`fuzzyand` <- function(...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    other.args <- list(...)
    
    negated <- FALSE
    
    if (is.vector(other.args[[1]])) {
        
        if (length(other.args) == 1) {
            return(other.args[[1]])
        }
        
        funargs <- toupper(unlist(lapply(match.call(expand.dots = TRUE), deparse)[-1]))
        negated <- grepl("1 -", funargs)
        other.args <- as.data.frame(other.args)
        tc <- capture.output(tryCatch(QCA::getName(funargs), error = function(e) e, warning = function(w) w))
        if (!grepl("simpleError", tc)) {
            names(other.args) <- QCA::getName(funargs)
        }
        else {
            names(other.args) <- LETTERS[seq(ncol(other.args))]
        }
    }
    else if (is.matrix(other.args[[1]]) | is.data.frame(other.args[[1]])) {
        other.args <- other.args[[1]]
    }
    
    cols <- colnames(other.args)
    
    
    if (is.null(cols)) {
        cols <- LETTERS[seq(ncol(other.args))]
    }
    
    if (length(negated) > 1) {
        cols[negated] <- tolower(cols[negated])
    }
    
    result <- apply(other.args, 1, min)
    attr(result, "name") <- paste(cols, collapse = "*")
    attr(result, "names") <- NULL
    return(structure(result, class = "fuzzyop"))
}
