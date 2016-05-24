`findSubsets` <-
function(noflevels, row.no, maximum) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    maximum <- ifelse(missing(maximum), prod(noflevels), maximum)
    result <- lapply(row.no, function(x) {
        # .Call("findSubsets", x, noflevels - 1, rev(c(1, cumprod(rev(noflevels))))[-1], maximum, PACKAGE="QCAGUI")
        QCA::callFindSubsets(x, noflevels - 1, rev(c(1, cumprod(rev(noflevels))))[-1], maximum)
    })
    return(sort(unique(unlist(result))))
}

