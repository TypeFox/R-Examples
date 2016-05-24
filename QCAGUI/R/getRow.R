`getRow` <- 
function(noflevels, row.no, zerobased=FALSE) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    max.combs <- prod(noflevels)
    if (any(row.no > (max.combs - zerobased))) {
        cat("\n")
        stop(simpleError(paste("There cannot be more than ", max.combs, " rows.\n\n", sep="")))
    }
    
    if (!zerobased) {row.no <- row.no - 1}
    
    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    t(sapply(row.no, function(x) x %/% mbase) %% noflevels)
}

