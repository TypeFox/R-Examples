summary.RSNPset.pvalue <- function(object, sort="p", decreasing=FALSE, nrows=10, dropcols=c(""), verbose=FALSE, ...) {
    if(nrows > attr(object, 'K')) {
        nrows <- attr(object, 'K')
    }
    cols <- setdiff(names(object),dropcols)
    if(sort %in% cols){
        rws <- order(object[,sort], decreasing = decreasing)[1:nrows]
    } else {
        rws <- order(row.names(object), decreasing = decreasing)[1:nrows]
    }
    
    if(verbose){
        cat("\n")
        if(attr(object, 'B') > 0){
            if(attr(object, 'pval.transform')) {
                cat("- Resampling p-values (pB, PB) come from comparison of\np-values across",attr(object, 'B'),"replications.\n")
            } else {
                cat("- Resampling p-values (pB) come from comparison of\ntest statistics across",attr(object, 'B'),"replications.\n")
            }
        }
        cat("- Q-values based on",attr(object, 'K'),"SNP sets.\n")
        cat('\n')
    }
    
    return(object[rws,cols])
}
