## appropriate behavior ???

## IF all missing data -- create multiPhylo4
## IF some have data -- create multiPhylo4d (user can coerce to multiPhylo4)
## IF (checkData) then stop if all data not identical to first data
##
## need constructors for multiPhylo4, multiPhylo4d!!
## FIXME: need code to construct tree.names ...

## function to bind trees together into a multi-tree object
tbind <- function(...,checkData=TRUE) {
    L <- list(...)
    namevec <- names(L)
    treeclasses <- c("multiPhylo4d","multiPhylo4","phylo4","phylo4d")
    tdataclasses <- c("multiPhylo4d","phylo4d")
    classes <- sapply(L,class)
    if (!all(classes %in% treeclasses)) {
        stop("all elements must be trees or multitrees")
    }
    hasData <- any(classes %in% tdataclasses)
    allData <- all(classes %in% tdataclasses)
    xfun <- function(x) {
        switch(class(x),
               phylo4=x,
               phylo4d=extractTree(x),
               multiPhylo4=x@phylolist,
               multiPhylo4d=suppressWarnings(as("multiPhylo4",x)@phylolist))}
    ## decompose multi-trees into lists
    treelist <- unlist(lapply(L,xfun))
    if (hasData) alldat <- lapply(L[classes %in% tdataclasses], tdata,
        type="tip")
    hasNodeData <- sapply(L[classes %in% tdataclasses], hasNodeData)
    if (any(hasNodeData)) warning("internal node data discarded")
    if (checkData) {
        ident <- sapply(alldat,identical,y=alldat[[1]])
        if (!all(ident)) stop(paste("tip data sets differ"))
    } ## ?? implement code to check which ones differ (taking
    ## null/multiple values in original set into account)
    if (hasData) return(new("multiPhylo4d",phylolist=treelist,
                            tip.data=alldat[[1]]))
    return(new("multiPhylo4",phylolist=treelist))
}
            

